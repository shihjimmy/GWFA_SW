import os
import argparse

BASE_BITS = 3  # 3 bits for base code: A/T/C/G/N -> 0..4

def truncate(file_path, NUM_EDGES):
    base_name = os.path.splitext(file_path)[0]
    out_path = f"{base_name}_truncate_to_{NUM_EDGES}.txt"

    # Map base to 3-bit code (0..4). Unknown -> 4 (N)
    base2code = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 'N': 4}

    # ----------------------------------------------------------------------
    # Parse GFA: collect node sequences and adjacency (out-edges)
    # ----------------------------------------------------------------------
    with open(file_path, "r") as f:
        lines = f.readlines()

    nodes = {}  # node_id -> sequence string
    edges = {}  # node_id -> list of neighbor node_ids

    for line in lines:
        parts = line.rstrip("\n").split('\t')
        if not parts:
            continue

        tag = parts[0]
        if tag == "S":
            # S   <id>    <sequence>   ...
            nid = int(parts[1])
            seq = parts[2]
            nodes[nid] = seq

        elif tag == "L":
            # L   <from>  <from_orient>  <to>  <to_orient>  <overlap>
            u = int(parts[1])
            v = int(parts[3])

            # We assume DAG (forward-only). If violated, we abort.
            # (Note: ID order != topological order; consider topo-sort in future.)
            if v < u:
                print("reverse edge encountered:", parts)
                raise RuntimeError("Graph is not forward-only/DAG as assumed.")

            edges.setdefault(u, []).append(v)

    # ----------------------------------------------------------------------
    # Build a global prefix position (flattened coordinate)
    # NOTE: This uses sorted node_ids. If GFA 'S' IDs are not already in
    # topological/linearized order, replace this with a topo-order.
    # ----------------------------------------------------------------------
    node_ids = sorted(nodes.keys())
    node_len = {nid: len(nodes[nid]) for nid in node_ids}

    prefix = {}  # node_id -> global start index (0-based) of its first base
    pos = 0
    for nid in node_ids:
        prefix[nid] = pos
        pos += node_len[nid]

    # ----------------------------------------------------------------------
    # Encode per-base lines:
    # - Top BASE_BITS: base code (A=0,T=1,C=2,G=3,N=4)
    # - Next NUM_EDGES bits (LSBs): distance mask, MSB=distance 1, LSB=distance NUM_EDGES
    #   For internal bases (within a node), always set distance-1 bit.
    #   For the last base of a node, set bits for each outgoing edge distance:
    #        d = prefix[v] - (prefix[u] + len(u) - 1)
    #   Only keep 1 <= d <= NUM_EDGES; larger distances are truncated.
    # ----------------------------------------------------------------------
    node_edge_bits = []

    for nid in node_ids:
        seq = nodes[nid]
        L = len(seq)

        # Internal bases: open distance-1 bit (stay within the same node)
        for i in range(L - 1):
            base = seq[i]
            base_code = base2code.get(base, 4)  # default to N (4) if unknown
            # put BASE_BITS at the top, distance mask at lower NUM_EDGES bits
            edge_bits = (base_code << NUM_EDGES) | (1 << (NUM_EDGES - 1))
            node_edge_bits.append(edge_bits)

        # Last base of the node: set outgoing distances
        base = seq[-1]
        base_code = base2code.get(base, 4)
        edge_bits = (base_code << NUM_EDGES)

        if nid in edges:
            for v in edges[nid]:
                # Distance from current node's last base to neighbor's first base
                d = prefix[v] - (prefix[nid] + L - 1)
                if 1 <= d <= NUM_EDGES:
                    edge_bits |= (1 << (NUM_EDGES - d))
                # If d > NUM_EDGES, it's not encodable in this scheme; ignore.

        node_edge_bits.append(edge_bits)

    # ----------------------------------------------------------------------
    # Write out: each line has BASE_BITS (base) + NUM_EDGES (distance mask) bits
    # ----------------------------------------------------------------------
    total_bits = NUM_EDGES + BASE_BITS
    with open(out_path, 'w') as out_file:
        for bits in node_edge_bits:
            line = bin(bits)[2:].zfill(total_bits)
            out_file.write(line + "\n")

    print(f"Edge bits have been written to {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert .gfa to truncated text (3-bit base + distance mask)")
    parser.add_argument("file", help="Path to input .gfa file")
    parser.add_argument("--num_edges", type=int, default=12,
                        help="Number of distance bits (default=12)")
    args = parser.parse_args()

    truncate(args.file, args.num_edges)


if __name__ == "__main__":
    main()
