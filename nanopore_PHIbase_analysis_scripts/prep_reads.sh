#!/bin/bash

#Script that takes raw reads and outputs concatenated fast reads for the barcode
#Also calcualtes contig stats for the barcode - output into different folder
#Filters the pass reads based on their length
#Run it within results directory -> ./prep_reads.sh barcode raw_read_location

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <barcode> <raw_read_location> <filter_length>"
    exit 1
fi

# Extract the arguments
barcode_number="$1"
location="$2"
filter_length="$3"

# Check if the specified location exists
if [ ! -d "$location" ]; then
    echo "Error: Location '$location' does not exist."
    exit 1
fi

# Create the barcode directory in current directory if it doesn't exist (-p flag)
barcode_dir="./barcode${barcode_number}"
mkdir -p "$barcode_dir"

# Concatenate the pass files based on the provided barcode number and location
zcat "$location/fastq_pass/barcode${barcode_number}"/* > "$barcode_dir/${barcode_number}_barcode.fastq"

#Check if the concatenated file exists
if [ ! -e "$barcode_dir/${barcode_number}_barcode.fastq" ]; then
    echo "Error: Concatenated FASTQ file '$barcode_dir/${barcode_number}_barcode.fastq' not found."
    exit 1
fi

# Run the Perl script (get_contig_stats) on the concatenated FASTQ file
get_contig_stats.pl -q -i "$barcode_dir/${barcode_number}_barcode.fastq" > $barcode_dir/${barcode_number}_contig_stats.txt

#Unzip fail barcodes 
zcat "$location/fastq_fail/barcode${barcode_number}"/* > "$barcode_dir/fail_barcode${barcode_number}.fastq"
wc -l $barcode_dir/fail_barcode${barcode_number}.fastq | tail -1 | awk '{print $1/4}' > $barcode_dir/${barcode_number}_num_fail.txt

# Delete the fail_barcode${barcode_number}.fastq file
rm "$barcode_dir/fail_barcode${barcode_number}.fastq"

#Filter pass reads on length and make new fastq file
awk -f /ei/projects/7/724b6a9a-6416-47eb-be58-d3737023e29b/scratch/getBigReads.awk -v min=${filter_length} $barcode_dir/${barcode_number}_barcode.fastq > $barcode_dir/${barcode_number}_barcode_${filter_length}bp.fastq

#Percentage of reads lost in filtering
echo "scale=2; (100 * $(wc -l < "$barcode_dir/${barcode_number}_barcode_${filter_length}bp.fastq") / $(wc -l < "$barcode_dir/${barcode_number}_barcode.fastq"))" | bc > "$barcode_dir/${barcode_number}_barcode_percent_retained.txt"
