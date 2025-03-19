#!/bin/bash
#SBATCH -J emerge_pipeline                    # Job name
#SBATCH -o logs/emerge_pipeline_sub_%a.out     # Output for each job
#SBATCH -e logs/emerge_pipeline_sub_%a.err     # Error for each job
#SBATCH -p ei-short                            # For testing on short
#SBATCH -c 1                                    # Number of CPU cores
#SBATCH --mem=2G                                # Memory size - optimise later
#SBATCH --mail-type=END,FAIL                    # Notifications for job done & fail 
#SBATCH --mail-user=mia.berelson@earlham.ac.uk  # Email address 
#SBATCH --array=01-96%20                        # Array range and concurrency limit (Set array number during submission)


#Usage------------
# num_barcodes=1
# sbatch array=1 submit_job_array.sh config_template.sh
#--------------

# Check if the configuration file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

# Source the configuration file
config=$1
if [ ! -f "$config" ]; then
    echo "Configuration file not found: $config"
    exit 1
fi
source "$config"

#Change directory to the output directory
cd $output_dir

# Create txt files if they don't exist
touch ./percent_reads_retained_length_filter.txt
touch ./no_fail_reads.txt
touch ./no_reads_ignored_parse_filter.txt
touch ./all_taxaID_count.tsv
touch ./lcaparse_summary.txt
touch ./lcaparse_perread.txt

# Add a header line to lcaparse_summary.txt if it doesn't already exist
grep -q "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank" ./lcaparse_summary.txt || printf "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank\n" > ./lcaparse_summary.txt
# Add a header line to lcaparse_perread.txt if it doesn't already exist
grep -q "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity" ./lcaparse_perread.txt || printf "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity\n" > ./lcaparse_perread.txt

# Get the barcode corresponding to this task ID
barcode_number="${barcode_list[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Now barcode_number holds the specific barcode to use for this task
echo "Processing barcode: $barcode_number"

# Set up directory for the current barcode
barcode_dir="barcode${barcode_number}"

# Execute the main processing script 
 /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/single_barcode_process.sh \
    $barcode_number $location $filter_length $reference_database $scratch_dir $output_dir $concatenated $contig_stats

# Check if the job script encountered an error and handle cancellation
if [ $? -ne 0 ]; then
    echo "Error detected for barcode ${barcode_number}. Cancelling job ${SLURM_JOB_ID}."
    scancel "$SLURM_JOB_ID"
fi
