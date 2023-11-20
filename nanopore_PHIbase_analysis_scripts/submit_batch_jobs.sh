#!/bin/bash

# Define an array of barcode numbers
barcode_list=("08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19")

# Define raw read location
location=(
    "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/raw/RL_24hCubTests_21092023/RL_24hCubTests_21092023/20230921_1707_X4_FAW72641_f261fc8c/")


#Define length to filter reads to 
filter_length=("300")

#Define mapping database
reference_database=(
    /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/phibase_rr_091123.fa)

# Create txt files if they don't exist
touch ./percent_reads_retained_length_filter.txt
touch ./no_fail_reads.txt
touch ./no_reads_ignored_parse_filter.txt
touch ./all_taxaID_count.tsv

# Submit SLURM jobs in a loop for each barcode
for barcode_number in "${barcode_list[@]}"; 
    do 
    
    # Calculate the barcode directory path
    barcode_dir="barcode${barcode_number}"
    #Create log direcrtory if it doesn't exist
    mkdir -p "$barcode_dir/logs"
    
    #Submit Jobs
    sbatch \
    --mem "2G" \
    -c 1 \
    -o "$barcode_dir/logs/${barcode_number}_submission_log.txt" \
    --error "$barcode_dir/logs/${barcode_number}_submission.err" \
    --job-name="${barcode_number}_submission" \
    --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/single_barcode_process.sh $barcode_number $location $filter_length $reference_database"
done
