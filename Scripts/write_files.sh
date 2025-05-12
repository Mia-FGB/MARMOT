#!/bin/bash
#script that adds different variables to files shared by the pipeline

# Check for the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <barcode>" 
    exit 1
fi

# Extract the arguments
barcode_number="$1"
barcode_dir="barcode${barcode_number}"

# Append percentage reads retained after filtering on length to file
echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_barcode_percent_retained.txt")" >> "percent_reads_retained_length_filter.txt"
echo "${barcode_number}_barcode_percent_retained.txt" added to "percent_reads_retained_length_filter.txt"

# Check if the number of failed reads file is empty and append appropriate message
if [ ! -s "${barcode_dir}/${barcode_number}_num_fail.txt" ]; then
    echo "Barcode_${barcode_number}: Fail reads not analysed" >> "no_fail_reads.txt"
else
    echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_num_fail.txt")" >> "no_fail_reads.txt"
fi

# Add a barcode column to the lcaparse_summary file
awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_summary.txt" >> "lcaparse_summary.txt"
echo "${barcode_number}_lcaparse_summary.txt" added to "lcaparse_summary.txt"

# Add a barcode column to the lcaparse_perread file
awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_perread.txt" >> "lcaparse_perread.txt"
echo "${barcode_number}_lcaparse_perread.txt" added to "lcaparse_perread.txt"

# Process the genome coverage file to add a barcode column and append only the first five columns (skipping readIDS)
input_file="${barcode_dir}/${barcode_number}_genome_coverage.txt"
output_file="genome_coverage_all.txt"

awk -v barcode="${barcode_number}" 'BEGIN {OFS="\t"} NR>1 {print barcode, $1, $2, $3, $4, $5}' "${input_file}" >> "${output_file}"
echo "${barcode_number}_genome_coverage.txt" added to "genome_coverage_all.txt"

# Combine the separate read numbers into a single file
cat "${barcode_dir}/${barcode_number}_read_no.tsv" >> read_numbers.tsv
echo "${barcode_number}_read_no.tsv" appended to "read_numbers.tsv"
