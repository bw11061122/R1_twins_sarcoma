#!/bin/bash

input_file="Data/chr1_388_reads_all.sam"

output_file=""
sample_name=""

while IFS= read -r line; do
    if [[ $line == Processing\ bam\ file\ * ]]; then
        sample_name=$(echo "$line" | awk -F'/' '{print $(NF-1)}')
        output_file="Data/chr1_388_reads_${sample_name}.txt"
        > "$output_file"
    else 
        echo "$line" >> "$output_file"
    fi
done < "$input_file"  