#!/usr/bin/env python3
import argparse
import time
import numpy as np

#######################################################################################
# Config
#######################################################################################

tile             = 256
continuous_match = 32
retreat_step     = 2

# band you want to use; they will be clamped to [-tile+1, tile-1]
set_lower        = -84
set_upper        =  43
lower            = min(0, set_lower)
upper            = max(0, set_upper)
band_lower       = max(-tile + 1, lower)
band_upper       = min(tile - 1, upper)

NUM_NODES  = tile
NUM_EDGES  = 12          # distance bits
BASE_BITS  = 3           # 3 bits for base (A/T/C/G/N)
code_to_base = {0: 'A', 1: 'T', 2: 'C', 3: 'G', 4: 'N', 5: ''}  # 5:'' sentinel

print(f"Tile size:          {tile}")
print(f"Retreat step:       {retreat_step}")
print(f"Band lower:         {band_lower}")
print(f"Band upper:         {band_upper}")
print(f"Continuous Match:   {continuous_match}")

#######################################################################################
# Helpers
#######################################################################################

def check_upper_lower(frontier):
    """
    Return (max_idx, min_idx, num_active) over the anti-diagonal frontier
    where entries are non-zero.
    """
    non_zero_indices = [i for i, v in enumerate(frontier) if v != 0]
    if not non_zero_indices:
        return None, None, 0
    return max(non_zero_indices), min(non_zero_indices), len(non_zero_indices)


def _seeds_from_frontier(frontier, band_lower, band_upper, NUM_NODES, nQ, nN):
    """
    Collect active (i, j) seeds from current frontier within band.
      frontier[k] = y+1 on anti-diagonal k = (i - j) + (NUM_NODES - 1)
    """
    seeds = []
    for i_off in range(band_lower, band_upper + 1):
        k = i_off + NUM_NODES - 1
        if 0 <= k < len(frontier) and frontier[k] != 0:
            j = frontier[k] - 1
            i = i_off + j
            if 0 <= i < nQ and 0 <= j < nN:
                seeds.append((i, j))
    return seeds


def _greedy_extend(i, j, query, nodes, edges,
                   frontier, traceback, stop_points,
                   NUM_NODES, NUM_EDGES, last,
                   branch_queue, out_border,
                   band_lower, band_upper):
    """
    Greedy diagonal extend using ONLY d=1 edges (pure ONE_MASK):
      - step (i+1, j+1) while eb == ONE_MASK and bases match and inside band
      - stop at boundary / mismatch / branching (eb != ONE_MASK)
      - if branching, enqueue candidate (i+1, j+d) within band that matches next base
      - if any branch mismatches next base, mark current (i,j) as stop_point
      - if eb == 0 (no outgoing), mark (i,j) as stop_point
      - if stepping beyond the last tile's nodes and not last tile, register out_border
    """
    ONE_MASK = 1 << (NUM_EDGES - 1)
    nQ, nN   = len(query), len(nodes)

    while True:
        # border of the tile
        if i == nQ - 1 or j == nN - 1:
            out_border.add((i, j))
            return

        eb = edges[j]
        if eb == 0:
            stop_points.add((i, j))
            return

        only_d1 = (eb == ONE_MASK)

        if only_d1:
            nj = j + 1
            if nj >= nN:
                if not last:
                    out_border.add((i, j))
                return

            # band check for (i+1, nj)
            delta = (i + 1) - nj
            if not (band_lower <= delta <= band_upper):
                # cannot continue inside band; current becomes a stop
                stop_points.add((i, j))
                return

            if query[i + 1] != nodes[nj]:
                stop_points.add((i, j))
                return

            # advance frontier & traceback at (i+1, nj)
            k = delta + NUM_NODES - 1
            if 0 <= k < len(frontier) and frontier[k] < nj + 1:
                frontier[k] = nj + 1
            if traceback[i + 1][nj] == "":
                traceback[i + 1][nj] = "1M"

            # keep walking along the diagonal
            i += 1
            j = nj
            continue

        # branching case: consider all d where bit is set
        exist_branch_fail = False

        for t in range(NUM_EDGES):  # MSB..LSB → d = 1..NUM_EDGES
            if not (eb & (1 << t)):
                continue

            d  = NUM_EDGES - t
            nj = j + d

            # stepping out of this tile's node slice
            if nj >= nN:
                if not last:
                    out_border.add((i, j))
                # do not return; check other branches
                continue

            delta = (i + 1) - nj
            if not (band_lower <= delta <= band_upper):
                continue

            if query[i + 1] == nodes[nj]:
                branch_queue.append((i + 1, nj, d))
            else:
                exist_branch_fail = True

        # after evaluating branches, stop here; extension will resume if a
        # branch is actually started (i.e., strictly advances the frontier)
        break

    if exist_branch_fail:
        # at least one candidate branch mismatched next base → this is a stop point
        stop_points.add((i, j))


def GWFA_retreat(traceback, x, y, edit_distance, step_config, continuous_match):
    """
    Retreat used before returning out-of-tile:
      - backtrack over consecutive 'M' (free matches), with continuous match cap
      - then retreat exactly one I/D/U, decreasing edit_distance
      - repeat up to step_config times or until hitting "0M"
    """
    retreat_step = step_config

    while retreat_step > 0:
        conti_match_count = 0
        old_x, old_y, old_edit = x, y, edit_distance

        # eat consecutive M
        while True:
            cell = traceback[x][y]
            last_move_pos = cell[:-1]
            last_move_dir = cell[-1]

            # start cell ("0M")
            if last_move_pos == "0" and last_move_dir == "M":
                return old_edit, (old_x, old_y)

            if last_move_dir != "M":
                break

            conti_match_count += 1
            if conti_match_count == continuous_match:
                return old_edit, (old_x, old_y)

            x -= 1
            y -= int(last_move_pos)

        # now one edit step back (I/D/U)
        cell = traceback[x][y]
        last_move_pos = cell[:-1]
        last_move_dir = cell[-1]
        edit_distance -= 1

        if last_move_dir == 'I':
            x -= 1
        elif last_move_dir == 'D':
            y -= int(last_move_pos)
        elif last_move_dir == 'U':
            x -= 1
            y -= int(last_move_pos)

        retreat_step -= 1

    return edit_distance, (x, y)



def GWFA_expand(traceback, frontier, stop_points, edges, nodes, query, last,
                NUM_NODES, NUM_EDGES, band_lower, band_upper):
    """
    Expand only from stop_points by one edit step (I/D/U), with band checks:
      - I: (x+1, y)
      - D: (x,   y+d)   for each outgoing edge bit
      - U: (x+1, y+d)
    Update frontier only if the new (y'+1) on that anti-diagonal is strictly farther
    than the current value. Write traceback the first time a cell is reached.
    """
    nR, nC = traceback.shape          # (rows, cols) = (NUM_NODES, NUM_NODES)
    L      = len(frontier)
    nN     = len(nodes)

    for (x, y) in stop_points:
        # I move: (x+1, y)
        if x + 1 < nR:
            delta_I = (x + 1) - y
            if band_lower <= delta_I <= band_upper:
                kI = delta_I + NUM_NODES - 1
                if 0 <= kI < L:
                    new_val = y + 1
                    if new_val > frontier[kI]:
                        frontier[kI] = new_val
                        if 0 <= y < nC and traceback[x + 1][y] == "":
                            traceback[x + 1][y] = "0I"

        # D/U moves from column y
        eb = edges[y] if 0 <= y < len(edges) else 0
        for t in range(NUM_EDGES):  # MSB..LSB → d = 1..NUM_EDGES
            if not (eb & (1 << t)):
                continue

            d      = NUM_EDGES - t
            next_y = y + d

            # last tile guard: do not step beyond local node slice
            if next_y >= nN and last:
                continue

            # if next character matches, greedy extension would have handled it.
            # Here we only materialize edits after a mismatch/break.
            if x + 1 < len(query) and 0 <= next_y < nN:
                if query[x + 1] == nodes[next_y]:
                    continue

            # D: (x, next_y)
            if 0 <= next_y < nC:
                delta_D = x - next_y
                if band_lower <= delta_D <= band_upper:
                    kD = delta_D + NUM_NODES - 1
                    if 0 <= kD < L:
                        new_val = next_y + 1
                        if new_val > frontier[kD]:
                            frontier[kD] = new_val
                            if traceback[x][next_y] == "":
                                traceback[x][next_y] = f"{d}D"

            # U: (x+1, next_y)
            if x + 1 < nR and 0 <= next_y < nC:
                delta_U = (x + 1) - next_y
                if band_lower <= delta_U <= band_upper:
                    kU = delta_U + NUM_NODES - 1
                    if 0 <= kU < L:
                        new_val = next_y + 1
                        if new_val > frontier[kU]:
                            frontier[kU] = new_val
                            if traceback[x + 1][next_y] == "":
                                traceback[x + 1][next_y] = f"{d}U"

#######################################################################################
# One-tile controller
#######################################################################################

def GWFA_tile_x_tile(nodes, edges, query, last,
                     NUM_NODES, NUM_EDGES, band_lower, band_upper,
                     continuous_match, retreat_step):
    """
    One tile controller:
      1) Extend all active anti-diagonals greedily along d=1, collect stop_points + candidate branches (band-limited).
      2) If any out-border → return (with retreat if not last).
      3) Iteratively start only branches that strictly advance the frontier; re-extend (band-limited).
      4) When no further branches can start → expand from stop_points (I/D/U), edit_distance += 1.
      5) Repeat.
    """
    nQ, nN = len(query), len(nodes)

    frontier  = [0] * (NUM_NODES + NUM_NODES - 1)
    traceback = np.full((NUM_NODES, NUM_NODES), "", dtype=object)

    # start cell
    frontier[NUM_NODES - 1] = 1
    traceback[0][0] = "0M"

    edit_distance = 0
    out_border    = set()

    while True:
        branch_queue = []
        stop_points  = set()

        # 1) greedy extend from current frontier
        seeds = _seeds_from_frontier(frontier, band_lower, band_upper,
                                     NUM_NODES, nQ, nN)
        for (xi, yj) in seeds:
            _greedy_extend(xi, yj, query, nodes, edges,
                           frontier, traceback, stop_points,
                           NUM_NODES, NUM_EDGES, last,
                           branch_queue, out_border,
                           band_lower, band_upper)

        # 2) out-border (tile exit)
        if out_border:
            upper, lower, num_active = check_upper_lower(frontier)
            # tie-break: prefer larger y; if y ties, prefer smaller x
            bx, by = max(out_border, key=lambda p: (p[1], p[0]))
            if not last:
                edit_distance, (bx, by) = GWFA_retreat(
                    traceback, bx, by, edit_distance, retreat_step, continuous_match
                )
            return edit_distance, traceback, (bx, by), (upper, lower, num_active)

        # 3) start branches that strictly improve frontier, then re-extend
        while True:
            started    = []
            seen_start = set()

            for (i1, j1, d) in branch_queue:
                k = i1 - j1 + NUM_NODES - 1
                if not (0 <= k < len(frontier)):
                    continue
                prev_off = frontier[k]        # old (y+1)
                if prev_off >= j1 + 1:
                    continue
                if (i1, j1) in seen_start:
                    continue

                # replace the frontier at this anti-diagonal with (j1+1)
                frontier[k] = j1 + 1
                if traceback[i1][j1] == "":
                    traceback[i1][j1] = f"{d}M"

                started.append((i1, j1))
                seen_start.add((i1, j1))

            if not started:
                break

            # re-extend from newly started seeds
            branch_queue = []
            for (xi, yj) in started:
                _greedy_extend(xi, yj, query, nodes, edges,
                               frontier, traceback, stop_points,
                               NUM_NODES, NUM_EDGES, last,
                               branch_queue, out_border,
                               band_lower, band_upper)

            if out_border:
                upper, lower, num_active = check_upper_lower(frontier)
                bx, by = max(out_border, key=lambda p: (p[1], p[0]))
                if not last:
                    edit_distance, (bx, by) = GWFA_retreat(
                        traceback, bx, by, edit_distance, retreat_step, continuous_match
                    )
                return edit_distance, traceback, (bx, by), (upper, lower, num_active)

        # 4) expand once from stop_points
        edit_distance += 1
        GWFA_expand(traceback, frontier, stop_points, edges, nodes, query, last,
                    NUM_NODES, NUM_EDGES, band_lower, band_upper)

#######################################################################################
# Tiled driver
#######################################################################################

def GWFA(fa_path, gfa_trunc_path):
    """
    Tiled driver:
      - Read truncated graph file (BASE_BITS + NUM_EDGES bits per line).
      - Read query FASTA (assumes sequence is on the second line).
      - Process tiles until either sequence or graph tile is exhausted.
      - Optionally force head-to-head / tail-to-tail end (to compare with original GWFA).
    """
    # graph (tile slice) container: sentinel first
    nodes = [code_to_base[5]]
    query = [code_to_base[5]]
    edges = [1 << (NUM_EDGES - 1)]

    # read truncated graph bitstrings: first BASE_BITS = base, then NUM_EDGES bits = mask
    with open(gfa_trunc_path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if len(s) != BASE_BITS + NUM_EDGES:
                raise ValueError(f"bitstring length {len(s)} != BASE_BITS+NUM_EDGES={BASE_BITS+NUM_EDGES}")
            
            base_code = int(s[:BASE_BITS], 2)
            nodes.append(code_to_base.get(base_code, 'N'))
            edges.append(int(s[BASE_BITS:], 2))

    
    # read query FASTA (sequence line)
    with open(fa_path, 'r') as f:
        lines = f.readlines()
        sequence = lines[1].strip()
        for b in sequence:
            query.append(b)

    print("-----------------------------")
    print(f"There are {len(nodes)} genome in the graph")
    print(f"There are {len(query)} genome in the sequence")
    print("-----------------------------")
    print("Start GWFA calculation.")
    print("-----------------------------")

    batch_size   = tile
    start_time   = time.time()
    x = y        = 0
    edit_total   = 0
    paths        = ""
    upper_max, lower_min, active_max = 0, tile * 2 - 1, 0

    while x < len(query) - 1 and y < len(nodes) - 1:
        path = ""

        batch_query = query[x : x + batch_size]
        batch_nodes = nodes[y : y + batch_size]
        batch_edges = edges[y : y + batch_size]

        left_x = len(query) - x - 1
        left_y = len(nodes) - y - 1
        last   = (left_x < batch_size and left_y < batch_size)

        score, tb, (end_x, end_y), (upper, lower, num_active) = GWFA_tile_x_tile(
            batch_nodes, batch_edges, batch_query, last,
            NUM_NODES, NUM_EDGES, band_lower, band_upper,
            continuous_match, retreat_step
        )

        # stuck guard
        if end_x == 0 and end_y == 0:
            if score != 0:
                print("ERROR: ended tile at (0,0) with non-zero score")
            break

        x += end_x
        y += end_y
        edit_total += score

        if upper is not None and upper > upper_max:
            upper_max = upper
        if lower is not None and lower < lower_min:
            lower_min = lower
        if num_active > active_max:
            active_max = num_active

        # reconstruct tile path (optional consistency check)
        # NOTE: here tb[end_x][end_y] stores only last step; if you need full
        # path consistency, keep your previous code.
        while True:
            if tb[end_x][end_y] == "0M":
                paths += path
                break
            step = tb[end_x][end_y]
            path = step + path
            last_move_pos = step[:-1]
            last_move_dir = step[-1]
            if last_move_dir == 'I':
                end_x -= 1
            elif last_move_dir == 'D':
                end_y -= int(last_move_pos)
            else:  # 'M' or 'U'
                end_x -= 1
                end_y -= int(last_move_pos)


    # (Optional) force head-to-head / tail-to-tail to compare with original GWFA
    if x != len(query)-1:
        diff_x = abs(len(query) - x - 1) 
        x = len(query) - 1
        edit_total += diff_x
        paths += diff_x * "0I"


    if y != len(nodes)-1:
        stack = [y]
        delete_finish = False
        diff_y = abs(len(nodes) - y -1)

        small_trace = np.full((diff_y,), "", dtype=object)
        
        while stack and not delete_finish:
            cur = stack.pop()
            edge = edges[cur]

            for i in range(NUM_EDGES):
                if edge & (1 << i):
                    next_cur = cur + (NUM_EDGES-i)

                    if next_cur < len(nodes)-1 and small_trace[next_cur - y - 1] == "":
                        stack.append(next_cur)
                        small_trace[next_cur - y - 1] = str(NUM_EDGES-i) + "D"
                        
                    elif next_cur == len(nodes)-1 and small_trace[next_cur - y - 1] == "":
                        delete_finish = True
                        small_trace[next_cur - y - 1] = str(NUM_EDGES-i) + "D"
                        y = len(nodes)-1
                        break
            
        
            if delete_finish:
                idx = diff_y-1
                delete_path = ""
                
                while True:
                    delete_path = small_trace[idx] + delete_path
                    edit_total += 1
                    
                    if idx == 0:
                        break

                    last_move_pos = small_trace[idx][0:-1]                  
                    idx -= int(last_move_pos)
                    
                    if idx < 0:
                        break
                    
                paths += delete_path



    elapsed = time.time() - start_time
    print(f"Execution time                               : {elapsed} seconds")
    print(f"Final Edit Distance                          : {edit_total}")
    print(f"Final Ending Position                        : {(x, y)}")
    print(f"wave upper band                              : {upper_max - NUM_NODES + 1}")
    print(f"wave lower band                              : {lower_min - NUM_NODES + 1}")
    print(f"Requiered wave band = upper-lower+1          : {upper_max - lower_min + 1}")
    print(f"Max active diagonal                          : {active_max}")
    print("-----------------------------")


    # quick consistency on total edits from path (optional)
    edit_count = sum(1 for ch in paths if ch in "UID")
    if edit_count == edit_total:
        print("traceback pass")
    else:
        print("traceback fail :", edit_count)
    print("-----------------------------")
    if x == len(query) - 1 and y == len(nodes) - 1:
        print("Ending Pass")
    else:
        print("Ending fail")


#######################################################################################
# Main
#######################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GWFA (tile x tile) on truncated graph")
    parser.add_argument('fa_file',             type=str, help="Path to the .fa file")
    parser.add_argument('truncated_gfa_file',  type=str, help="Path to the truncated gfa file (.txt)")
    args = parser.parse_args()
    GWFA(args.fa_file, args.truncated_gfa_file)
