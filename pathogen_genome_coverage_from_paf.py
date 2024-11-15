#!/usr/bin/python

#A script to calculate genome coverage for specific genera from minimap output
#Filters on identity and coverage
#Should be run from directory which contains barcode directories

import sys, getopt, errno, os, csv

#This block of code means it can take any paf file 
# When using the -b flag in the terminal will take barcodes
try:
    opts, args = getopt.getopt(sys.argv[1:],"b:")
except getopt.GetoptError:
    print("Option not recognised.")
    print("python my_script.py -b <barcode_number>")
    sys.exit(2)
for opt, arg in opts:
    if opt == "-b":
        input_arg = arg
    else:
        print("python my_script.py -i <barcode_number>")
        sys.exit(2)
barcode_number = input_arg

def main(barcode_number):
    pafFilename = None  # Define pafFilename with a default value
    barcode_dir = None 

    # Set the directory where output files should be saved
    barcode_dir = "./barcode{}".format(barcode_number)


    # Combine the barcode directory and the PAF file name to create the full path
    pafFilename = os.path.join(barcode_dir, "{}_mapped.paf".format(barcode_number))

    # Check if pafFilename is an empty string before attempting to open the file
    if not pafFilename:
        print("pafFilename is not valid. Exiting.")
        sys.exit(2)

main(barcode_number)

def calculate_coverage(pafFilename):
    read_longest_alignment = {} # Dictionary to store the longest alignment for each read
    ignored_reads = [] # list to store readIDs of ignored reads

    with open(pafFilename, 'r') as paf_file:
        reader = csv.reader(paf_file, delimiter='\t')

        for row in reader:
            # Extract relevant fields
            read_id = row[0]                # Name of the read (column 1)
            q_length = int(row[1])           # Length of the read (column 2)
            q_start = int(row[2])            # Start position of the alignment in the read (column 3)
            q_end = int(row[3])              # End position of the alignment in the read (column 4)
            taxaID =  row[5].split("|")[1] #for PHIbase 
            #taxaID =  row[5] #for other fasta references
            mq = int(row[11])                # Mapping quality (column 12)
            matching_bases = int(row[9])     # Number of matching bases (column 10)
            total_bases = int(row[10])       # Length -> Number bases, including gaps, in the mapping (col 11)
            identity = (matching_bases / total_bases)* 100 # Calculate the identity of the mapping
            coverage =((q_end - q_start) / q_length) * 100 # Calculate the coverage of the mapping

        # Filter reads on coverage and identity
            if total_bases >=100 and identity >= 70: #still working on this
                # Check if the read is already encountered
                if read_id not in read_longest_alignment:
                    # If new, add the read and its matching base count
                    read_longest_alignment[read_id] = matching_bases
                else:
                    # If already seen, store the highest number of matching bases
                    if matching_bases > read_longest_alignment[read_id]:
                        read_longest_alignment[read_id] = matching_bases

            else: #not quite right as this doesn't consider reads with multiple alignments that have been lost in loop above
                # Add the read_id to the ignored_reads list and more info 
                ignored_reads.append((read_id, identity, coverage))

    
    total_matching_bases = sum(read_longest_alignment.values())
    number_of_mapped_reads = len(read_longest_alignment) #to work out how many unique read IDs were considered 'mapped' to the reference