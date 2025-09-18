#!/bin/bash

# Clear out all files inside directories starting with out_trim_
for dir in ./out_trim_err*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        echo "Clearing contents of directory: $dir"
        rm -rf "$dir"/*  # Remove all files inside the directory but keep the directory
    fi
done


for dir in ./pbsim3_trim_err*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        echo "Clearing contents of directory: $dir"
        rm -rf "$dir"/*  # Remove all files inside the directory but keep the directory
    fi
done


# Iterate through all subdirectories starting with "out_sequence_"
for dir in ./out_sequence_err*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        len=$(basename "$dir" | sed 's/out_sequence_//')  # Get the part after "out_sequence_"
        echo "Processing directory: $dir for sequence with length: $len"
        python3 GWFA_gen_testdata.py 1  $len
        python3 GWFA_gen_testdata.py 22 $len
    fi
done

python3 GWFA_change_all_edge.py

echo "All test cases have been processed."
