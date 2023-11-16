#!/bin/bash

for dir in $(find . -mindepth 1 -type d); do
    files_1=("$dir"/*_1.fq.gz)
    if [ $(echo "${#files_1[@]}") -gt 1 ]; then
        echo "Concatenating ${#files_1[@]} files in $dir"
        cat "${files_1[@]}" > "$dir/$(basename "${files_1[0]}" | awk -F'_' '{print $1"_"$2}')_concatenated_1.fq.gz"
        echo "Concatenation complete for FORWARD reads in $dir"
    fi

    files_2=("$dir"/*_2.fq.gz)
    if [ $(echo "${#files_2[@]}") -gt 1 ]; then
        echo "Concatenating ${#files_2[@]} files in $dir"
        cat "${files_2[@]}" > "$dir/$(basename "${files_2[0]}" | awk -F'_' '{print $1"_"$2}')_concatenated_2.fq.gz"
        echo "Concatenation complete for REVERSE reads in $dir"
    fi
done

