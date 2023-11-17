#!/bin/bash

echo -e "Sample_ID\tsR1\tsR2" > samples.tsv

find /hdd0/raid3/scratch/MEGshared/ea_biofilm/data -type f \( -name "*_1.fq.gz" -o -name "*_2.fq.gz" \) | sort | \
while read -r file; do
    dir=$(dirname "$file")
    sample_id=$(basename "$dir")

    if [[ $file == *_1.fq.gz ]]; then
        sR1="$file"
    elif [[ $file == *_2.fq.gz ]]; then
        sR2="$file"
        echo -e "$sample_id\t$sR1\t$sR2" >> samples.tsv
    fi
done


# Remove duplicate lines based on sample ID and add sample_ID column
awk -F'\t' '!seen[$1]++ { OFS="\t"; if ($3 ~ /concatenated_[12]\.fq\.gz$/) print $1, $2, $3; else if (!($1 in samples)) { samples[$1] = 1; print $1, $2, $3 } }' samples.tsv > samples_no_duplicates.tsv

