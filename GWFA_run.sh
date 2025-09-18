#!/bin/bash

# 總計時開始
total_start=$SECONDS

# Iterate over all subdirectories starting with "out_trim_err"
for dir in ./out_trim_err*; do
    if [ -d "$dir" ]; then
        echo "Processing data in directory: $dir"
        dir_start=$SECONDS   # 子目錄計時開始

        len=$(basename "$dir" | sed 's/out_trim_//')
        TEST_CASE_FOLDER="$dir"
        OUTPUT_FILE="./GWFA_result_all/GWFA_${len}.txt"  

        # 如果結果檔已存在就刪掉
        [ -f "$OUTPUT_FILE" ] && rm "$OUTPUT_FILE"

        # 遍歷目錄下所有 .fa 檔案
        for fa_file in "$TEST_CASE_FOLDER"/*.fa; do
            base_name=$(basename "$fa_file" .fa)
            gfa_file="$TEST_CASE_FOLDER/${base_name}_truncate_to_12.txt"

            if [ -f "$gfa_file" ]; then
                echo "Running test case pair: $fa_file and $gfa_file"
                echo "Running test case pair: $fa_file and $gfa_file" >> "$OUTPUT_FILE"

                python3 GWFA.py "$fa_file" "$gfa_file" >> "$OUTPUT_FILE" 2>&1

                # 分隔線
                printf '%0.s-' {1..50} >> "$OUTPUT_FILE"
                echo >> "$OUTPUT_FILE"
                printf '%0.s ' {1..50} >> "$OUTPUT_FILE"
                echo >> "$OUTPUT_FILE"
                printf '%0.s-' {1..50} >> "$OUTPUT_FILE"
                echo >> "$OUTPUT_FILE"
            else
                echo "Warning: No matching .gfa file found for $fa_file"
            fi
        done

        # 子目錄計時結束
        dir_end=$SECONDS
        dir_duration=$((dir_end - dir_start))
        dir_minutes=$((dir_duration / 60))
        dir_seconds=$((dir_duration % 60))

        # 將時間資訊輸出到終端與檔案
        echo "All test cases for $dir processed. Results saved to $OUTPUT_FILE."
        echo "Time spent on $dir: ${dir_minutes} minutes ${dir_seconds} seconds."
        echo "----------------------------------------------"

        echo >> "$OUTPUT_FILE"
        echo "Time spent on $dir: ${dir_minutes} minutes ${dir_seconds} seconds." >> "$OUTPUT_FILE"
        echo "----------------------------------------------" >> "$OUTPUT_FILE"
    fi
done

# 總計時結束
total_end=$SECONDS
total_duration=$((total_end - total_start))
total_minutes=$((total_duration / 60))
total_seconds=$((total_duration % 60))

echo "All test cases have been processed."
echo "Total execution time: ${total_minutes} minutes ${total_seconds} seconds."

