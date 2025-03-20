# Configuration file for the pipeline

# Pipeline parameters
location="/ei/projects/7/78b740c8-8ef2-451f-a465-48c035f07ada/data/raw/rebasecalled_reads/RL_MaizeFirstWeek_24022025_sup"
# location="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/raw/RL_24hCubTests_21092023/RL_24hCubTests_21092023/20230921_1707_X4_FAW72641_f261fc8c"

filter_length="300"

reference_database="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/pathogen_database_080524.fa"

# Directory for temporary files
scratch_dir="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/pipeline_test"

# Directory for output files
output_dir="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/pipeline_test"

# List of barcodes to process
# Uncomment the following line to process a sequence of barcodes
# barcode_list=($(seq -w 01 88))
# Uncomment and modify the following line to process only specific barcodes
barcode_list=("11" "02" "05" "08")

# Set to "yes" if the reads are already concatenated into one file, else set to "no"
concatenated="yes"

# Set to "yes" if you want to calculate contig stats for the barcode, else set to "no"
contig_stats="yes"

genome_lengths_file="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/080524_genome_lengths.tsv"