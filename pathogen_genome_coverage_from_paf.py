#!/usr/bin/python

#A script to calculate genome coverage for specific genera from minimap output
#Filters on identity and coverage
#Should be run from directory which contains barcode directories
#Uses a genome_lengths_file generated with genome_lengths_from_fasta.py script using the reference database used for mapping

import sys, getopt, errno, os, csv

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python pathogen_genome_coverage_from_paf.py <barcode_number> <genome_lengths_file>")
    sys.exit(1)

# Extract arguments
barcode_number = sys.argv[1]
genome_lengths_file = sys.argv[2]

# Set the directory where output files should be saved
barcode_dir = "./barcode{}".format(barcode_number)

# Combine the barcode directory and the PAF file name to create the full path
pafFilename = os.path.join(barcode_dir, "{}_mapped.paf".format(barcode_number))

# Check if pafFilename is an empty string before attempting to open the file
if not pafFilename:
    print("pafFilename is not valid. Exiting.")
    sys.exit(2)

genome_coverage_file = os.path.join(barcode_dir, "{}_genome_coverage.txt".format(barcode_number))
excluded_reads_file = os.path.join(barcode_dir, "{}_coverage_excluded_reads.txt".format(barcode_number))
multi_taxa_reads_file = os.path.join(barcode_dir, "{}_coverage_multi_taxa_reads.txt".format(barcode_number))

def calculate_coverage(pafFilename, genome_lengths_file, genome_coverage_file, ignored_reads_file):
    # Dictionary to store total mapped bases for each taxaID
    taxa_mapped_bases = {}
    # List to store information about ignored reads
    excluded_reads = []
    # Dictionary to track reads that map to multiple taxa
    multi_taxa_reads = {}
    # Set to track processed read_id and taxaID pairs
    processed_reads = set()

    # Load the genome lengths into a dictionary
    genome_lengths = {}  # Use a separate variable to store genome lengths
    with open(genome_lengths_file, 'r') as genome_lengths_f:
        reader = csv.reader(genome_lengths_f, delimiter='\t')
        next(reader)  # Skip header row
        for row in reader:
            taxaID = row[0]
            length = int(row[1])  # Genome length in base pairs
            genome_lengths[taxaID] = length

    with open(pafFilename, 'r') as paf_file:
        reader = csv.reader(paf_file, delimiter='\t')

        for row in reader:
            # Extract relevant fieldss
            read_id = row[0]                 # Name of the read (column 1)
            q_length = int(row[1])           # Length of the whole read (column 2)
            a_start = int(row[2])            # Start position of the alignment in the read (column 3)
            a_end = int(row[3])              # End position of the alignment in the read (column 4)
            taxaID =  row[5].split("|")[1]   #for PHIbase references
            mq = int(row[11])                # Mapping quality (column 12)
            matching_bases = int(row[9])     # Number of matching bases (column 10)
            a_length = int(row[10])          # Length -> Number bases, including gaps, in the mapping (col 11)
            identity = (matching_bases / a_length)* 100 # Calculate identity of mapping = how many bases map / length of alignment
            coverage =((a_length) / q_length) * 100 # Calculate coverage of mapping = length of alignment compared to read
           
            # Filter reads on coverage and identity thresholds
            if coverage >= 80 and identity >= 70:  # Can change these values
                # Check if this read_id and taxaID pair has already been processed
                if (read_id, taxaID) not in processed_reads:
                    # Add the q_length to the taxa_mapped_bases for this taxaID
                    if taxaID not in taxa_mapped_bases:
                        taxa_mapped_bases[taxaID] = 0
                    taxa_mapped_bases[taxaID] += q_length
                    processed_reads.add((read_id, taxaID))

                    # Track reads that map to multiple taxa and pass filters
                    if read_id not in multi_taxa_reads:
                        multi_taxa_reads[read_id] = []
                    multi_taxa_reads[read_id].append((taxaID, identity, coverage))

            else:
                # Track alignments which do not pass the filters
                excluded_reads.append((read_id, taxaID, q_length, a_length, identity, coverage))
    
    # Calculate coverage using the mapped_bases and genome_lengths
    with open(genome_coverage_file, "w") as mapped_file:
        mapped_file.write("taxaID\tmapped_bases\tgenome_length\tcoverage_percentage\n")
        for taxaID, bases in taxa_mapped_bases.items():
            genome_length = genome_lengths.get(taxaID, 0)  # Default to 0 if taxaID not found
            if genome_length > 0:
                coverage_percentage = (bases / genome_length) * 100
            else:
                print(f"Warning: taxaID {taxaID} not found in genome lengths table.")
                coverage_percentage = 'N/A'
            mapped_file.write(f"{taxaID}\t{bases}\t{genome_length}\t{coverage_percentage:.4f}\n")

    # Write excluded reads to file
    with open(excluded_reads_file, "w") as excluded_file:
        excluded_file.write("read_id\ttaxaID\tread_length\talignment_length\tidentity\tcoverage\n")
        for read_id, taxaID, q_length, a_length, identity, coverage in excluded_reads:
            excluded_file.write(f"{read_id}\t{taxaID}\t{q_length}\t{a_length}\t{identity:.2f}\t{coverage:.2f}\n")

    # Write reads that map to multiple taxa to file
    with open(multi_taxa_reads_file, "w") as multi_taxa_file:
        multi_taxa_file.write("read_id\ttaxaIDs\tidentity\tcoverage\n")
        for read_id, taxaIDs in multi_taxa_reads.items():
            if len(taxaIDs) > 1:  # Only include reads that map to multiple taxa
                for taxaID, identity, coverage in taxaIDs:
                    multi_taxa_file.write(f"{read_id}\t{taxaID}\t{identity:.2f}\t{coverage:.2f}\n")
    
    print("Processed genome coverage for barcode")

calculate_coverage(pafFilename, genome_lengths_file, genome_coverage_file, excluded_reads_file)