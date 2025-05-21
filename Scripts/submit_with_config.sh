#!/bin/bash

# Usage: ./submit_with_config.sh config.sh

# Check for config
if [ $# -ne 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

config=$1
echo "Loading config file: $1"
source "$config"

# Check required variables exist
if [ -z "$sample" ] ; then
  echo "Your config file must define 'sample'"
  exit 1
fi

if [ -z "${barcode_list[*]}" ]; then
  echo "Your config file must define a 'barcode_list' array"
  exit 1
fi

# Compute the array range from the number of barcodes
num_barcodes=${#barcode_list[@]}
array_range=$(printf "01-%02d" "$num_barcodes")
echo "Detected $num_barcodes barcodes"
echo "Submitting SLURM array range: $array_range"

# Create logs directory if it doesn't exist
log_dir="logs/${sample}"
mkdir -p "$log_dir"
echo "Log files will be written to: $log_dir"

# Create a temporary SLURM script
temp_script=$(mktemp)
echo "Generating temporary SLURM script: $temp_script"


cat > "$temp_script" <<EOF
#!/bin/bash
#SBATCH -J $sample
#SBATCH -o $log_dir/submission_%a.out
#SBATCH -e $log_dir/submission_%a.err
#SBATCH -p ei-short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mia.berelson@earlham.ac.uk
#SBATCH --array=$array_range%96

echo "Running SLURM task ID: \$SLURM_ARRAY_TASK_ID"
echo "Loaded job for sample: $sample"
echo "Using config file: $(realpath "$config")"

# Run the actual job script
bash $(which submit_job_array.sh) $(realpath "$config")
EOF

# Show the SLURM script - uncomment for debugging
# echo "Final SLURM submission script:"
# echo "-------------------------------"
# cat "$temp_script"
# echo "-------------------------------"

# Submit the job and capture the job ID
array_job_id=$(sbatch "$temp_script" | awk '{print $4}')
echo "Submitted array job ID: $array_job_id"

# Once the array job is submitted, we can submit the write_files job
# This job will write the shared output files after all array jobs are done
write_files_job_id=$(sbatch --dependency=afterok:$array_job_id \
    --mem 1G \
    -p ei-short \
    -o "$log_dir/write_files.out" \
    -e "$log_dir/write_files.err" \
    --job-name="${sample}_write_files" \
    --wrap "bash $(which write_files.sh) $(realpath "$config")" | awk '{print $4}')
echo "Submitted write files job ID: $write_files_job_id"


# Create the risk plots job, run after write_files_job has completed
create_plots_id=$(sbatch --dependency=afterok:$write_files_job_id  \
    --mem=5G \
    -p ei-short \
    -o "$log_dir/riskplots.out" \
    --error "$log_dir/riskplots.err" \
    --job-name="riskplots" \
    --wrap "source activate r-marmot_env && Rscript $(realpath create_risk_plots.R) $output_dir $risk_table_file" | awk '{print $4}')


echo "Submitted create plots job: $create_plots_id, will run when write_files_job ($write_files_job_id) finishes"

