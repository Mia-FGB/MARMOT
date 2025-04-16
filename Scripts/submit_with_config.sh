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
#SBATCH --array=$array_range%20

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

# Submit the job
sbatch "$temp_script"
