#!/bin/bash
# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <barcode_number> <location> <filter_length> <reference_database>"
    exit 1
fi

# Extract the arguments - these will be provided by the submit_batch_jobs.sh script
barcode_number="$1"
location="$2"
filter_length="$3"
reference_database="$4"

# Check if the specified location exists
if [ ! -d "$location" ]; then
    echo "Error: Location '$location' does not exist."
    exit 1
fi

# Calculate the barcode directory path
barcode_dir="barcode${barcode_number}"
#Create log directory if it doesn't exist
mkdir -p "$barcode_dir/logs"

#Name & create the files for script
percent_results_file="percent_reads_retained_length_filter.txt"
touch -c "$percent_results_file"

fail_number_file="no_fail_reads.txt"
touch -c "$fail_number_file"

ignored_reads_number_file="no_reads_ignored_parse_filter.txt"
touch -c "$ignored_reads_number_file"

all_taxaID_count="all_taxaID_count.tsv"
touch -c "$all_taxaID_count"

# Execute the prep_reads.sh script with the barcode and location
JOBID1=$(sbatch --mem 30G \
    -p ei-medium \
    -o "barcode${barcode_number}/logs/${barcode_number}_prepreads.out" \
    --error "barcode${barcode_number}/logs/${barcode_number}_prepreads.err" \
    --job-name="${barcode_number}_prepreads" \
    --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/prep_reads.sh $barcode_number $location $filter_length" | awk '{print $NF}')

echo "Submitted prep_reads job ($JOBID1)"

#Minimap Job after JOB1 is finished
JOBID2=$(sbatch --dependency=afterok:$JOBID1 \
    --mem 60G \
    -p ei-medium \
    -o "$barcode_dir/logs/${barcode_number}_minimap.out" \
    --error "$barcode_dir/logs/${barcode_number}_minimap.err" \
    --job-name="${barcode_number}_minimap2" \
    --wrap "minimap2 -x map-ont ${reference_database} \"$barcode_dir/${barcode_number}_barcode_${filter_length}bp.fastq\" \\
    > \"$barcode_dir/${barcode_number}_mapped.paf\"" | awk '{print $NF}')

echo "Submitted minimap job ($JOBID2) will run when ($JOBID1) finishes"
echo "Barcode ${barcode_number} is being mapped to ${reference_database}"

#Submit job for paf_parse script
JOBID3=$(sbatch \
    --dependency=afterok:$JOBID2 \
    --mem 30G \
    -p ei-medium \
    -o "$barcode_dir/logs/${barcode_number}_paf_parse.out" \
    --error "$barcode_dir/logs/${barcode_number}_paf_parse.err" \
    --job-name="${barcode_number}_paf_parse" \
    --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/paf_parse.py -b ${barcode_number}" | awk '{print $NF}')

echo "Submitted paf_parse ($JOBID3) will run when ($JOBID2) finishes"
echo "Barcode ${barcode_number} paf file is being parsed"

# Submit a new job with dependency on JOBID3
JOBID4=$(sbatch --dependency=afterok:$JOBID3 \
    --mem 2G \
    -p ei-medium \
    -o "$barcode_dir/logs/${barcode_number}_write_files_output.txt" \
    -e "$barcode_dir/logs/${barcode_number}_write_files_error.txt" \
    --job-name="${barcode_number}_write_files" \
    --wrap '
    # Append the barcode and percentage retained to the results file
    echo "Barcode_${barcode_number}: $(cat "barcode${barcode_number}/${barcode_number}_barcode_percent_retained.txt")" >> "$percent_results_file"

    # Print the result
    echo "Barcode ${barcode_number}: Percentage retained appended to $percent_results_file"

    # Append the barcode and number of fail reads to file
    echo "Barcode_${barcode_number}: $(cat "barcode${barcode_number}/${barcode_number}_num_fail.txt")" >> "$fail_number_file"

    # Print the result
    echo "Barcode ${barcode_number}: number of fail reads appended to $fail_number_file"

    # Append the barcode and number of ignored reads to file
    echo "barcode${barcode_number}: $(cat "barcode${barcode_number}/${barcode_number}_number_ignored_reads.tsv")" >> "$ignored_reads_number_file"

    # Print the result
    echo "Barcode ${barcode_number}: Number ignored reads appended to $ignored_reads_number_file"

    # Append the taxaID_count file to all taxaID_count file
    echo "$(cat "barcode${barcode_number}/${barcode_number}_taxaID_counts.tsv")" >> "$all_taxaID_count"

    # Print the result
    echo "Barcode ${barcode_number}: "${barcode_dir}/${barcode_number}_taxaID_counts.tsv" Taxa IDs and counts appended to $all_taxaID_count"
    ')

echo "Submitted write_files job with JOBID4: $JOBID4"


