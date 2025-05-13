#!/bin/bash

#Script that takes raw reads and outputs concatenated fast reads for the barcode
#Also calculates contig stats for the barcode - output into different folder
#Filters the pass reads based on their length
#Run it within results directory -> ./prep_reads.sh barcode raw_read_location

# Check for the correct number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage prep_reads.sh: $0 <barcode> <raw_read_location> <filter_length> <scratch_dir> <output_dir> <concatenated> <contig_stats>" 
    exit 1
fi

# Extract the arguments
barcode_number="$1"
location="$2"
filter_length="$3"
scratch_dir="$4"
output_dir="$5"
concatenated="$6"
contig_stats="$7"

# Check if the specified location exists
if [ ! -d "$location" ]; then
    echo "Error: Location '$location' does not exist."
    exit 1
fi

# Specify barcode directory - created in current directory 
barcode_dir="$output_dir/barcode${barcode_number}"

#Check if the flag file exists and exit without error if it does
if [ -f "$barcode_dir/prep_reads_finished.flag" ]; then
    echo "Flag file prep_reads_finished.flag already exists for barcode ${barcode_number}. Exiting"
    exit 0
fi

# Create scratch directory if it doesn't exist
mkdir -p "$scratch_dir"

#Only for non-concatenated reads
if [ "$concatenated" != "yes" ]; then

    # Concatenate the pass files based on the provided barcode number and location 
    # only unzips when necessary and outputs the new concatenated file into scratch area
    # First looks for fastq_pass directory 
    if [ -d "$location/fastq_pass/barcode${barcode_number}" ]; then
        for file in "$location/fastq_pass/barcode${barcode_number}"/*; do
            if [[ "$file" == *.fastq || "$file" == *.fq ]]; then               #Check if file already unzipped
                cat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
            elif [[ "$file" == *.gz ]]; then                                    #Unzip if needed
                zcat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
            else
                echo "Warning: Skipping unrecognized file format '$file'."
            fi
        done
        # Same loop as above but looking in fastq directory
    elif [ -d "$location/fastq/barcode${barcode_number}" ]; then
        for file in "$location/fastq/barcode${barcode_number}"/*; do
            if [[ "$file" == *.fastq || "$file" == *.fq ]]; then 
                cat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
            elif [[ "$file" == *.gz ]]; then 
                zcat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
            else
                echo "Warning: Skipping unrecognized file format '$file'."
            fi
        done
    else
        echo "Error: Neither fastq or fastq_pass directory exists for location ${location} and barcode ${barcode_number}."
        exit 1
    fi

    # Check if the fastq_fail directory exists
    if [ -d "$location/fastq_fail/barcode${barcode_number}" ]; then
        # Unzip fail barcodes
        zcat "$location/fastq_fail/barcode${barcode_number}"/* > "$scratch_dir/fail_barcode${barcode_number}.fastq"
        if [ ! -s "$scratch_dir/fail_barcode${barcode_number}.fastq" ]; then
            echo "Warning: No fail barcodes found for barcode ${barcode_number}."
        else
            # Count the number of fail barcodes
            if ! wc -l "$scratch_dir/fail_barcode${barcode_number}.fastq" | tail -1 | awk '{print $1/4}' > "$barcode_dir/${barcode_number}_num_fail.txt"; then
                echo "Error: Failed to count the number of fail barcodes."
            fi
        fi
    else
        echo "Warning: fastq_fail directory for barcode ${barcode_number} does not exist. Skipping fail barcode processing."
    fi

    # Filter pass reads on length and create a new fastq file
    if ! awk -f /ei/projects/7/724b6a9a-6416-47eb-be58-d3737023e29b/scratch/getBigReads.awk -v min="${filter_length}" "$scratch_dir/${barcode_number}_barcode.fastq" > "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq"; then
        echo "Error: Failed to filter pass reads based on length."
        exit 1
    fi

    # Calculate the percentage of reads retained after filtering
    if ! echo "scale=2; (100 * $(wc -l < "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq") / $(wc -l < "$scratch_dir/${barcode_number}_barcode.fastq"))" | bc > "$barcode_dir/${barcode_number}_barcode_percent_retained.txt"; then
        echo "Error: Failed to calculate the percentage of reads retained after filtering."
        exit 1
    fi

    #Process contig stats if requested on the concatenated file
    if [ "$contig_stats" == "yes" ]; then
        echo "Calculating contig stats for barcode ${barcode_number}..."
        # Run the Perl script (get_contig_stats) on the concatenated FASTQ file - not length filtered
        # Outputs the contig stats into the barcode directory not scratch area
        if ! get_contig_stats.pl -q -i "$scratch_dir/${barcode_number}_barcode.fastq" > "$barcode_dir/${barcode_number}_contig_stats.txt"; then
            echo "Error: Failed to run get_contig_stats.pl on the concatenated FASTQ file."
            exit 1
        fi
    fi
fi  # Ensures this block is skipped when concatenated = 'yes'

if [ "$concatenated" == "yes" ]; then
    # Check if the barcode directory exists
    if [ ! -d "$location/fastq/barcode${barcode_number}" ]; then
        echo "Error: Barcode directory '$location/fastq/barcode${barcode_number}' does not exist."
        exit 1
    fi

    # Check there is only one .fastq file in the barcode directory
    num_files=$(find "$location/fastq/barcode${barcode_number}" -maxdepth 1 -type f -name "*.fastq" | wc -l)
    if [ "$num_files" -ne 1 ]; then
        echo "Error: More than one .fastq file found in the barcode directory, not concatenated. Exiting without processing."
        exit 1
    fi

    # Get the actual filename
    file=$(find "$location/fastq/barcode${barcode_number}" -maxdepth 1 -type f -name "*.fastq")

    # Filter these concatenated reads on length & create new fasta file in scratch
    if ! awk -f /ei/projects/7/724b6a9a-6416-47eb-be58-d3737023e29b/scratch/getBigReads.awk -v min="${filter_length}" "$file" > "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq"; then
        echo "Error: Failed to filter pass reads based on length."
        exit 1
    fi

    # Calculate the percentage of reads retained after filtering
    total_reads=$(($(wc -l < "$file") / 4))
    retained_reads=$(($(wc -l < "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq") / 4))

    if ! echo "scale=2; (100 * $retained_reads / $total_reads)" | bc > "$barcode_dir/${barcode_number}_barcode_percent_retained.txt"; then
        echo "Error: Failed to calculate the percentage of reads retained after filtering."
        exit 1
    fi

    #Process contig stats if requested on the pre-existing file
    if [ "$contig_stats" == "yes" ]; then
        echo "Calculating contig stats for barcode ${barcode_number}..."
        # Run the Perl script (get_contig_stats) on the concatenated FASTQ file - not length filtered
        # Outputs the contig stats into the barcode directory not scratch area
        if ! get_contig_stats.pl -q -i "$file" > "$barcode_dir/${barcode_number}_contig_stats.txt"; then
            echo "Error: Failed to run get_contig_stats.pl on the concatenated FASTQ file."
            exit 1
        fi
    fi
fi


if [ "$contig_stats" != "yes" ]; then
    echo "Not calculating contig stats for original file of barcode ${barcode_number}..."
fi

# Run contig stats on the filtered reads, whether or not requested
echo "Calculating contig stats for $filter_length filtered barcode ${barcode_number}..."
# Run the Perl script (get_contig_stats) on the filtered FASTQ file
# Outputs the contig stats into the barcode directory not scratch area
if ! get_contig_stats.pl -q -i "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq" > "$barcode_dir/${barcode_number}_contig_stats_filtered_${filter_length}.txt"; then
    echo "Error: Failed to run get_contig_stats.pl on the filtered FASTQ file."
    exit 1
fi

#Write a tsv file with the number of reads before and after filtering
output_tsv="$barcode_dir/${barcode_number}_read_no.tsv"
touch "$output_tsv"
# Append the current barcode's data
echo -e "${barcode_number}\t${total_reads}\t${retained_reads}" >> "$output_tsv"

# Create the flag file to indicate successful completion
echo "prep_reads.sh has finished successfully, writing prep_reads_finished.flag"
touch "$barcode_dir/prep_reads_finished.flag"
