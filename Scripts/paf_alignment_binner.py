#!/usr/bin/python

# A script to take the PAF alignment file
# Look at the alignment distribution across the genome per taxa
# Bin the alignments 

import sys, os, csv


def process_paf_file(pafFilename, bin_count): #Function to process paf file one row at a time
    p_ram_binned = {}
    h_frax_binned = {}

    seen_reads_p_ram = set()
    seen_reads_h_frax = set()

    with open(pafFilename, 'r') as paf_file:
        reader = csv.reader(paf_file, delimiter='\t')
        for row in reader:
            read_id = row[0]
            q_length = int(row[1])
            ref_start = int(row[7])
            ref_end = int(row[8])
            taxaID = row[5].split("|")[1]
            matching_bases = int(row[9])
            a_length = int(row[10])
            identity = (matching_bases / a_length) * 100
            #coverage = ((a_length) / q_length) * 100

            if taxaID == "164328" and a_length >= 150 and identity >= 85: 
                # Phytophthora ramorum
                midpoint = (ref_start + ref_end) // 2 
                genome_length = 57451392
                p_ram_bin_size = genome_length // bin_count
                bin_index = (midpoint // p_ram_bin_size) * p_ram_bin_size
    
                key = (read_id, bin_index)
                # Use a set to track unique reads for Phytophthora ramorum
                # to avoid counting the same read multiple times in the same bin
                if key not in seen_reads_p_ram:
                    seen_reads_p_ram.add(key)
                    if bin_index not in p_ram_binned:
                        p_ram_binned[bin_index] = 0
                    p_ram_binned[bin_index] += 1

            elif taxaID == "746836" and a_length >= 150 and identity >= 85: 
                # Hymenoscyphus fraxineus
                midpoint = (ref_start + ref_end) // 2 
                genome_length = 63496644
                h_frax_bin_size = genome_length // bin_count
                bin_index = (midpoint // h_frax_bin_size) * h_frax_bin_size
    
                key = (read_id, bin_index)
                if key not in seen_reads_h_frax:
                    seen_reads_h_frax.add(key)
                    if bin_index not in h_frax_binned:
                        h_frax_binned[bin_index] = 0
                    h_frax_binned[bin_index] += 1

    return p_ram_binned, p_ram_bin_size, h_frax_binned, h_frax_bin_size

def write_binned_output(filename, binned_data, taxon_label, bin_size):
    with open(filename, "w") as out:
        out.write(f"taxon\tbin_start\tbin_end\tread_count\n")
        for bin_start in sorted(binned_data):
            count = binned_data[bin_start]
            bin_end = bin_start + bin_size
            out.write(f"{taxon_label}\t{bin_start}\t{bin_end}\t{count}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: paf_alignment_binner.py <barcode_number> <output_dir> <bin_count>")
        sys.exit(1)

    barcode_number = sys.argv[1]
    output_dir = sys.argv[2]
    bin_count = int(sys.argv[3])
    barcode_dir = os.path.join(output_dir, f"barcode{barcode_number}")

    print (f"Processing barcode {barcode_number} in directory {barcode_dir}")
    
    pafFilename = os.path.join(barcode_dir, f"{barcode_number}_mapped.paf")
    if not os.path.isfile(pafFilename):
        print(f"PAF file not found: {pafFilename}")
        sys.exit(2)

    print(f"Processing genome coverage for barcode{barcode_number} from {pafFilename}")

    # Process and bin alignments
    p_ram_binned, p_ram_bin_size, h_frax_binned, h_frax_bin_size = process_paf_file(pafFilename, bin_count)

    # Set output file paths
    p_ram_output = os.path.join(barcode_dir, f"{barcode_number}_p_ram_binned_{bin_count}.txt")
    h_frax_output = os.path.join(barcode_dir, f"{barcode_number}_h_frax_binned_{bin_count}.txt")

    # Write results
    write_binned_output(p_ram_output, p_ram_binned, "Phytophthora_ramorum", p_ram_bin_size)
    write_binned_output(h_frax_output, h_frax_binned, "Hymenoscyphus_fraxineus", h_frax_bin_size)

    print("Done: binned alignment summaries written.")

if __name__ == "__main__":
    main()