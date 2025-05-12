#!/bin/bash
# Script that adds different variables to files shared by the pipeline
# Runs just once at the end of the pipeline, each barcode appended sequentially to avoid overlap 

# Check if the configuration file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

config="$1"
if [ ! -f "$config" ]; then
    echo "Configuration file not found: $config"
    exit 1
fi

# Source the config file to load barcode_list
source "$config"

if [ -z "${barcode_list[*]}" ]; then
    echo "No barcode list found in config."
    exit 1
fi

if [ -z "$output_dir" ]; then
    echo "Output directory not specified in config."
    exit 1
fi

cd $output_dir

# Create txt files if they don't exist
touch ./percent_reads_retained_length_filter.txt
touch ./no_fail_reads.txt
touch ./lcaparse_summary.txt
touch ./lcaparse_perread.txt
touch ./genome_coverage_all.txt
touch ./read_numbers.tsv

# Add a header line to lcaparse_summary.txt if it doesn't already exist
grep -q "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank" ./lcaparse_summary.txt || printf "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank\n" > ./lcaparse_summary.txt
# Add a header line to lcaparse_perread.txt if it doesn't already exist
grep -q "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity" ./lcaparse_perread.txt || printf "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity\n" > ./lcaparse_perread.txt
# Add a header line to genome_coverage_all.txt if it doesn't already exist
grep -q "Barcode\ttaxaID\tmapped_bases\tgenome_length\tcoverage_percentage\tnum_reads" ./genome_coverage_all.txt || printf "Barcode\ttaxaID\tmapped_bases\tgenome_length\tcoverage_percentage\tnum_reads\n" > ./genome_coverage_all.txt
# Add a header line to read_numbers.tsv if it doesn't already exist
grep -q "Barcode\tRead_Count\tFilterReadCount" ./read_numbers.tsv || printf "Barcode\tRead_Count\tFilterReadCount\n" > ./read_numbers.tsv



# Loop through each barcode
for barcode_number in "${barcode_list[@]}"; do
    echo "Processing barcode: $barcode_number"
    barcode_dir="barcode${barcode_number}"

    echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_barcode_percent_retained.txt")" >> "percent_reads_retained_length_filter.txt"
    echo "${barcode_number}_barcode_percent_retained.txt" added to "percent_reads_retained_length_filter.txt"

    if [ ! -s "${barcode_dir}/${barcode_number}_num_fail.txt" ]; then
        echo "Barcode_${barcode_number}: Fail reads not analysed" >> "no_fail_reads.txt"
    else
        echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_num_fail.txt")" >> "no_fail_reads.txt"
    fi

    awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_summary.txt" >> "lcaparse_summary.txt"
    echo "${barcode_number}_lcaparse_summary.txt" added to "lcaparse_summary.txt"

    awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_perread.txt" >> "lcaparse_perread.txt"
    echo "${barcode_number}_lcaparse_perread.txt" added to "lcaparse_perread.txt"

    input_file="${barcode_dir}/${barcode_number}_genome_coverage.txt"
    output_file="genome_coverage_all.txt"
    awk -v barcode="${barcode_number}" 'BEGIN {OFS="\t"} NR>1 {print barcode, $1, $2, $3, $4, $5}' "${input_file}" >> "${output_file}"
    echo "${barcode_number}_genome_coverage.txt" added to "genome_coverage_all.txt"

    cat "${barcode_dir}/${barcode_number}_read_no.tsv" >> read_numbers.tsv
    echo "${barcode_number}_read_no.tsv" appended to "read_numbers.tsv"
done

