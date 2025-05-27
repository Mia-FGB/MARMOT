#!/bin/bash

# Script to process raw reads for a given barcode
# version 2 - code has been simplified using chatGPT from the original prep_reads.sh script
# This script concatenates reads, filters them by length, and calculates contig stats.
# Usage: ./prep_reads.sh <barcode> <input_dir> <min_length> <scratch_dir> <output_dir> <is_concatenated> <run_stats>

if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <barcode> <input_dir> <min_length> <scratch_dir> <output_dir> <is_concatenated> <run_stats>"
    exit 1
fi

barcode="$1"
input_dir="$2"
min_length="$3"
scratch="$4"
output="$5"
is_concat="$6"
run_stats="$7"

barcode_dir="$output/barcode${barcode}"
mkdir -p "$scratch" "$barcode_dir"

if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' not found."
    exit 1
fi

concat_fastq="$scratch/${barcode}_reads.fastq"
filtered_fastq="$scratch/${barcode}_reads_${min_length}bp.fastq"

if [ "$is_concat" != "yes" ]; then
    pass_dir=""
    for sub in fastq_pass fastq; do
        if [ -d "$input_dir/$sub/barcode${barcode}" ]; then
            pass_dir="$input_dir/$sub/barcode${barcode}"
            break
        fi
    done

    if [ -z "$pass_dir" ]; then
        echo "Error: No valid fastq_pass or fastq directory for barcode $barcode"
        exit 1
    fi

    for file in "$pass_dir"/*; do
        if [[ "$file" =~ \.fastq$ || "$file" =~ \.fq$ ]]; then
            cat "$file" >> "$concat_fastq"
        elif [[ "$file" == *.gz ]]; then
            zcat "$file" >> "$concat_fastq"
        else
            echo "Skipping unrecognised file: $file"
        fi
    done

    fail_dir="$input_dir/fastq_fail/barcode${barcode}"
    if [ -d "$fail_dir" ]; then
        zcat "$fail_dir"/* > "$scratch/fail_${barcode}.fastq"
        if [ -s "$scratch/fail_${barcode}.fastq" ]; then
            wc -l "$scratch/fail_${barcode}.fastq" | awk '{print $1/4}' > "$barcode_dir/${barcode}_num_fail.txt"
        else
            echo "No fail reads found for barcode $barcode"
        fi
    else
        echo "No fastq_fail directory for barcode $barcode"
    fi
else
    concat_fastq=$(find "$input_dir/fastq/barcode${barcode}" -maxdepth 1 -type f -name "*.fastq" ! -name "*.fastq.*")
    if [ $(echo "$concat_fastq" | wc -l) -ne 1 ]; then
        echo "Error: Expected exactly one .fastq file in barcode dir."
        exit 1
    fi
fi

awk -f /ei/projects/7/724b6a9a-6416-47eb-be58-d3737023e29b/scratch/getBigReads.awk -v min="$min_length" "$concat_fastq" > "$filtered_fastq"

total_reads=$(($(wc -l < "$concat_fastq") / 4))
retained_reads=$(($(wc -l < "$filtered_fastq") / 4))
echo "scale=2; (100 * $retained_reads / $total_reads)" | bc > "$barcode_dir/${barcode}_percent_retained.txt"

if [ "$run_stats" == "yes" ]; then
    get_contig_stats.pl -q -i "$concat_fastq" > "$barcode_dir/${barcode}_contig_stats.txt"
fi

echo "Running contig stats for filtered reads..."
get_contig_stats.pl -q -i "$filtered_fastq" > "$barcode_dir/${barcode}_contig_stats_filtered_${min_length}.txt"

echo -e "${barcode}\t${total_reads}\t${retained_reads}" >> "$barcode_dir/${barcode}_read_no.tsv"
echo "Script finished for barcode ${barcode}."
