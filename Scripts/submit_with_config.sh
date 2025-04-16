#!/bin/bash

# Usage: ./submit_with_config.sh config.sh

# Check for config
if [ $# -ne 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

config=$1
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

# Create logs directory if it doesn't exist
mkdir -p logs/${job_name}

# Create a temporary SLURM script
temp_script=$(mktemp)

cat > "$temp_script" <<EOF
#!/bin/bash
#SBATCH -J $job_name
#SBATCH -o logs/${job_name}/submission_%a.out
#SBATCH -e logs/${job_name}/submission_%a.err
#SBATCH -p ei-short
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mia.berelson@earlham.ac.uk
#SBATCH --array=$array_range%20

bash $(realpath submit_job_array.sh) $config
EOF

# Submit the generated script
sbatch "$temp_script"
