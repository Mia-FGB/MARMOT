# Configuration file for MARMOT Pipeline

# Pipeline parameters
sample="marmot_cf_2023"

# Location is the path to raw / rebasecalled reads 
# Should be up to the level above fastq directory
location="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/church_farm_2023/SubSamp_200k_CF_2023_Wkly_Comb_Reads_SUP"

# Filter length - reads shorter than this will eb removed from analysis 
filter_length="300"

# Path to the reference database for alignment 
reference_database="/ei/projects/8/818a2389-cfbb-4ffb-b424-98250a8a3118/data/results/PHIbase_database_May25/pathogen_database_042025.fa"

# Directory for temporary files
scratch_dir="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/MARMOT_Output/cf_2023"

# Directory for output files
output_dir="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/MARMOT_Output/cf_2023"

# List of barcodes to process
# Uncomment the following line to process a sequence of barcodes
barcode_list=($(seq -w 01 39))
# Uncomment and modify the following line to process only specific barcodes
# barcode_list=("11" "02" "05" "08")

# Set to "yes" if the reads are already concatenated into one file (will be the case if they have been rebasecalled), else set to "no"
concatenated="yes"

# Set to "yes" if you want to calculate contig stats for the barcode, else set to "no"
contig_stats="no"

# File that contains taxaIDS and genome lenghts, is created at the same time as the reference database
genome_lengths_file="/ei/projects/8/818a2389-cfbb-4ffb-b424-98250a8a3118/data/results/PHIbase_database_May25/042025_genome_lengths.tsv"

#Risk Table File - generated with generate_risk_table.py at the same time as the reference database
risk_table_file="/ei/projects/8/818a2389-cfbb-4ffb-b424-98250a8a3118/data/results/PHIbase_database_May25/risk_table.csv"

# Barcode Labels for the graphs as a tab-separated file
barcode_labels="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/MARMOT_Output/cf_2023/barcode_labels.tsv"