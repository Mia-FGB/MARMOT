#!/usr/bin/python

#A script to calculate genome coverage for specific genera from minimap output
#Filters on identity and coverage
#Should be run from directory which contains barcode directories

import sys, getopt, errno, os, csv

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python pathogen_genome_coverage_from_paf.py <barcode_number> <lengths_table>")
    sys.exit(1)

# Extract arguments
barcode_number = sys.argv[1]
lengths_table = sys.argv[2]

# Set the directory where output files should be saved
barcode_dir = "./barcode{}".format(barcode_number)

# Combine the barcode directory and the PAF file name to create the full path
pafFilename = os.path.join(barcode_dir, "{}_mapped.paf".format(barcode_number))

# Check if pafFilename is an empty string before attempting to open the file
if not pafFilename:
    print("pafFilename is not valid. Exiting.")
    sys.exit(2)

genome_coverage_file = os.path.join(barcode_dir, "{}_genome_coverage.txt".format(barcode_number))
ignored_reads_file = os.path.join(barcode_dir, "{}_ignored_reads.txt".format(barcode_number))

def calculate_coverage(pafFilename, genome_coverage_file, ignored_reads_file):
    # Nested dictionary to track longest alignment for each read_id per taxaID
    read_longest_alignment = {}
    # Dictionary to store total mapped bases for each taxaID
    taxa_mapped_bases = {}
    # List to store information about ignored reads
    ignored_reads = []

    with open(pafFilename, 'r') as paf_file:
        reader = csv.reader(paf_file, delimiter='\t')

        for row in reader:
            # Extract relevant fieldss
            read_id = row[0]                # Name of the read (column 1)
            q_length = int(row[1])           # Length of the read (column 2)
            q_start = int(row[2])            # Start position of the alignment in the read (column 3)
            q_end = int(row[3])              # End position of the alignment in the read (column 4)
            taxaID =  row[5].split("|")[1] #for PHIbase references
            mq = int(row[11])                # Mapping quality (column 12)
            matching_bases = int(row[9])     # Number of matching bases (column 10)
            total_bases = int(row[10])       # Length -> Number bases, including gaps, in the mapping (col 11)
            identity = (matching_bases / total_bases)* 100 # Calculate the identity of the mapping
            coverage =((q_end - q_start) / q_length) * 100 # Calculate the coverage of the mapping

            # Filter reads on coverage and identity thresholds
            if coverage >= 80 and identity >= 70: #can change these values
                # Initialize read_id entry in read_longest_alignment if not already present
                if read_id not in read_longest_alignment:
                    read_longest_alignment[read_id] = {}
                
                # Check if this taxaID has been encountered for this read_id
                if taxaID not in read_longest_alignment[read_id]:
                    # If not present, initialize with current alignment length
                    read_longest_alignment[read_id][taxaID] = total_bases
                else:
                    # Update to keep the longest alignment for this taxaID
                    if total_bases > read_longest_alignment[read_id][taxaID]:
                        read_longest_alignment[read_id][taxaID] = total_bases

            else:
                # Add the read_id to the ignored_reads list with additional information
                ignored_reads.append((read_id, taxaID, total_bases, identity, coverage))

    # Aggregate the total alignment lengths per taxaID
    for read_id, taxa_dict in read_longest_alignment.items():
        for taxaID, longest_alignment in taxa_dict.items():
            if taxaID not in taxa_mapped_bases:
                taxa_mapped_bases[taxaID] = 0
            taxa_mapped_bases[taxaID] += longest_alignment

    
    #Write mapped bases per taxaID to file
    with open(genome_coverage_file, "w") as mapped_file:
        mapped_file.write("taxaID\tmapped_bases\n")
        for taxaID, bases in taxa_mapped_bases.items():
            mapped_file.write(f"{taxaID}\t{bases}\n")

    # Write ignored reads information to file
    with open(ignored_reads_file, "w") as ignored_file:
        ignored_file.write("read_id\ttaxaID\ttotal_bases\tidentity\tcoverage\n")
        for read_id, taxaID, total_bases, identity, coverage in ignored_reads:
            ignored_file.write(f"{read_id}\t{taxaID}\t{total_bases}\t{identity:.2f}\t{coverage:.2f}\n")

    print(f"Results written to {genome_coverage_file} and {ignored_reads_file}")

calculate_coverage(pafFilename, genome_coverage_file, ignored_reads_file)