#!/bin/bash

#Usage------------
# num_barcodes=1
# sbatch array=1 submit_job_array.sh config_template.sh
# This script is combined with a wrapper to submit a job array
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

echo "In submit_job_array.sh"
echo "Config file path: $config"
echo "Sourcing config file..."

source "$config"
echo "Loaded config for sample: $sample"

# Make output directory if it doesn't exist
mkdir -p "$output_dir"
#Change directory to the output directory
cd $output_dir
echo "Data will be written to: $output_dir"

# Get the barcode corresponding to this task ID
barcode_number="${barcode_list[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Now barcode_number holds the specific barcode to use for this task
echo "Processing barcode: $barcode_number"

# Set up directory for the current barcode
barcode_dir="barcode${barcode_number}"

echo "Running single_barcode_process.sh with:"
echo "$barcode_number $location $filter_length $reference_database $scratch_dir $output_dir $concatenated $contig_stats $genome_lengths_file"

# Execute the main processing script 
 /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/single_barcode_process.sh \
    "$barcode_number" "$location" "$filter_length" "$reference_database" "$scratch_dir" "$output_dir" "$concatenated" "$contig_stats" "$genome_lengths_file"

# Check if the job script encountered an error and handle cancellation
if [ $? -ne 0 ]; then
    echo "Error detected for barcode ${barcode_number}. Cancelling job ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}."
    scancel "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    echo "Job for barcode ${barcode_number} completed successfully."
fi