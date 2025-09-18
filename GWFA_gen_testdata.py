import argparse
import random
from tqdm import tqdm
import re


#############################################################################

""" 
    Change this for more generated testdata
"""
num_of_data = 1000 # for each chromosome
code_to_base = {0: 'A', 1: 'T', 2: 'C', 3: 'G', 4: ' '}
#############################################################################


def remove_dashes(sequence):
    # Remove all dashes from the sequence
    return sequence.replace("-", "")


def binary_search(arr, target):
    left, right = 0, len(arr) - 1

    while left+1 < right:
        mid = (left + right) // 2
        # print(left, right, mid)
        
        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid
        else:
            right = mid

    return left



def read_and_split_file(file_path, delimiter=" "):
    result = []
    
    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split(delimiter)
            result.append(split_line)
    
    return result


def extract_gfa_data(gfa_lines, start, pos):
    nodes = []
    edges = []
    
    for line in gfa_lines:
        if line.startswith('S'):  
            parts = line.split()
            node_id = int(parts[1])
            
            if start <= node_id <= pos:
                nodes.append(parts)
                
        elif line.startswith('L'): 
            parts = line.split()
            node1 = int(parts[1])
            node2 = int(parts[3])

            if node1 >= start and node2 <= pos:
                edges.append(parts)

    nodes.sort(key=lambda x: int(x[1])) 
    edges.sort(key=lambda x: (int(x[1]), int(x[3])))

    return nodes, edges



def write_to_file(start_pos_data, gfa_lines, chrom, num_of_data, length):
    random.seed(120)
    sampled_indices = random.sample(range(len(start_pos_data)), min(num_of_data, len(start_pos_data)))

    with open(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_trim.txt", "r") as f:
        fq_lines = f.readlines()


    # Create a file to save sampled indices
    sampled_indices_file = f"./out_trim_{length}/chr{chrom}_sampled_indices.txt"
    

    # Write the sampled indices to the file
    with open(sampled_indices_file, 'w') as idx_file:
        for idx in sampled_indices:
            idx_file.write(f"{idx+1}\n")
    
    # Confirm that sampled indices have been saved
    print(f"Sampled indices saved to {sampled_indices_file}")
    

    for idx in sampled_indices:

        start_pos = start_pos_data[idx] 
        name = start_pos[0]

        start = int(start_pos[1])
        pos = int(start_pos[3])
        start_offset = int(start_pos[2])
        end_offset = int(start_pos[4])

        nodes_in_range, edges_in_range = extract_gfa_data(gfa_lines, start, pos)

        file_name = f"./out_trim_{length}/chr{chrom}_{name}.gfa"
        sequence_file_name = f"./out_trim_{length}/chr{chrom}_{name}.fa" 

        with open(file_name, 'w') as f:
            
            """
                Get the subgraph and cutting off the part we don't need      
            """

            for i in range(len(nodes_in_range)):
                
                if i==0:
                    nodes_in_range[i][2] = nodes_in_range[i][2][start_offset:]

                elif i==len(nodes_in_range)-1:
                    nodes_in_range[i][2] = nodes_in_range[i][2][:end_offset]

                f.write("\t".join(str(x) for x in nodes_in_range[i]) + "\n")
                    
            
            for edge in edges_in_range:
                f.write("\t".join(str(x) for x in edge) + "\n")
        

        with open(sequence_file_name, 'w') as seq_file:
            sequence = fq_lines[4*idx + 1].strip()
            seq_file.write(f">{name}\n{sequence}\n")
        

        print(f"Data for start={name} written to {file_name} and sequence written to {sequence_file_name}")

print("---------------------------------------------")



#"---------------------------------------------""---------------------------------------------""---------------------------------------------"


parser = argparse.ArgumentParser(description="generate_gfa_fa_truncated_gfa")
parser.add_argument('chrom', type=str, help="chromosome")
parser.add_argument('len', type=str, help="length of pbsim3 sequence")

args     = parser.parse_args()
chrom    = args.chrom
length   = args.len


match = re.search(r'(\d+)(k|K)', length)  # This will match the number and the 'k'/'K'

if match:
    # Extract the numeric part and convert to int
    numeric_value = int(match.group(1))
    
    # Format with 'K' if it's not already
    formatted_length = f"{numeric_value}k"
    #print(formatted_length)
else:
    print("Invalid length format")


f = open(f"./out_sequence_{length}/chr{chrom}_pbsim3_{formatted_length}_0001.fq", "r")
f2 = open(f"./out_sequence_{length}/chr{chrom}_pbsim3_{formatted_length}_0001.maf", "r")
lines = f.readlines()
lines_maf = f2.readlines()

f.close()
f2.close()

f = open(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_trim.txt", "w", newline='\n') 
f2 = open(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_start_pos.txt", "w", newline='\n')

num_reads = len(lines_maf)//4
after_trim = 1

for i in tqdm(range(num_reads), desc="Processing"):

    if "N" not in lines[4*i+1]:
        lines_maf[(4*i)+1] = [i for i in lines_maf[(4*i)+1].split(" ") if i != '']
        lines_maf[(4*i)+2] = [i for i in lines_maf[(4*i)+2].split(" ") if i != '']

        # remove character beyond "ATCG"
        lines_maf[(4*i)+2][6] = "".join([c for c in lines_maf[(4*i)+2][6] if c in "ATCG\n"])
        lines[(4*i)+3] = lines[(4*i)+3][:len(lines[(4*i)+1])-1] + "\n"
        
        f.writelines([f"@S1_{after_trim}\n", lines_maf[(4*i)+2][6], f"+S1_{after_trim}\n", lines[(4*i)+3]])
        f2.write(f"S1_{after_trim} {lines_maf[(4*i)+1][2]} {lines_maf[(4*i)+1][3]} {lines_maf[(4*i)+2][4]}\n")
        after_trim += 1


f.close()
f2.close()
print("fq_trim is finished")

print("---------------------------------------------")

# read graphs
f = open(f"./out_graph/chr{chrom}.gfa", "r")
lines = f.readlines()
f.close()

nodes = dict()

for i in tqdm(lines[1::], desc="Processing Lines"):
    line = i[0:-1].split()

    if line[0] == "S":
        nodes[int(line[1])] = line[2]
        
    # path of main sequence
    elif line[0] == "P":
        path = line[2].split("+,")
        path[-1] = path[-1][:-1]



# position of node on reference
path_accu_len = [0 for i in range(len(path))]
temp = 0

for i in range(len(path)):
    path_accu_len[i] = temp
    temp += len(nodes[int(path[i])])


f = open(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_start_pos.txt", "r")
lines = f.readlines()
f.close()

f = open(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_pos_on_graph.txt", "w")

for i in tqdm(lines, desc="Processing"):
    i = i.split(" ")

    start_node = binary_search(path_accu_len, int(i[1]))
    start_offset = int(i[1]) - path_accu_len[start_node]

    end_node = binary_search(path_accu_len, int(i[1])+int(i[2]))
    end_offset = int(i[1])+int(i[2]) - path_accu_len[end_node]

    if end_offset == 0:
        # Perfectly fits the node before the already found end_node
        end_node -= 1
        end_offset = path_accu_len[end_node] - path_accu_len[end_node-1]
    
    f.write(f"{i[0]} {path[start_node]} {start_offset} {path[end_node]} {end_offset} {i[3]}")


f.close()
print("fq_pos is finished")

print("---------------------------------------------")

f = open(f"./out_graph/chr{chrom}.gfa", "r")
gfa_lines = f.readlines()
f.close()

start_pos_data = read_and_split_file(f"./pbsim3_trim_{length}/pbsim3_chr{chrom}_pos_on_graph.txt")
write_to_file(start_pos_data, gfa_lines, chrom, num_of_data, length)
print("gfa_trim is finished")









# print("---------------------------------------------")


# f = open(f"./out_graph/chr{chrom}.gfa", "r")
# lines = f.readlines()
# f.close()

# nodes = dict()
# edges = dict()


# for i in tqdm(lines, desc="Processing"):
#     line = i[0:-1].split('\t')

#     if line[0] == "S":
#         nodes[int(line[1])] = line[2]

#     elif line[0] == "L":
#         start = int(line[1])
#         end   = int(line[3])

#         # after vg's generation, it should be a DAG
#         if (end < start):
#             print("reverse edge", line)
#             break
#         try:
#             edges[start].append(end)
#         except:
#             edges[start] = [end]

# print("Finish reading in the whole genome graph.")
# print("Start to do the truncation.")

# for node in edges:
#     edges[node].sort()
    
    
# edge_vectors = {}

# for node, neighbors in edges.items():
#     edge_vector = 0  

#     for idx, neighbor in enumerate(neighbors):
#         diff = abs(neighbor - node)

#         if diff <= NUM_EDGES:
#             edge_vector |= (1 << (NUM_EDGES-diff))

#     edge_vectors[node] = edge_vector


# sorted_nodes = sorted(nodes.items()) 


# out_file = f"./out_graph/chr{chrom}_truncated.gfa"

# with open(out_file, "w") as out_file:
#     for node, value in sorted_nodes:
#         edge_vector = edge_vectors.get(node, 0) 
#         out_file.write(f"{node} {value} {bin(edge_vector)[2:].zfill(NUM_EDGES)}\n")  
        
        
# print("truncate the whole genome graph is finished.")
# print("---------------------------------------------")
