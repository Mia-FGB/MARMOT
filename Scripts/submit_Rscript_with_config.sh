#!/bin/bash

# Usage: ./submit_Rscript_with_config.sh config.sh

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

if [ -z "$output_dir" ]; then
  echo "Your config file must define 'output_dir'"
  exit 1
fi
if [ -z "$risk_table_file" ]; then
  echo "Your config file must define 'risk_table_file'"
  exit 1
fi

if [ -z "${barcode_labels[*]}" ]; then
  echo "Your config file must contain the path to the barcode labels file in 'barcode_labels'"
  exit 1
fi

# Create logs directory if it doesn't exist
log_dir="logs/${sample}"
mkdir -p "$log_dir"
echo "Log files will be written to: $log_dir"

# Need the full path to the R script here 
create_plots_id=$(sbatch \
    --mem=5G \
    -p ei-short \
    -o "$log_dir/riskplots_rep.out" \
    --error "$log_dir/riskplots_rep.err" \
    --job-name="${sample}_riskplots" \
    --wrap "source activate r-marmot_env && Rscript /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/create_risk_plots.R $output_dir $risk_table_file $barcode_labels" | awk '{print $4}')


echo "Submitted create plots job: $create_plots_id using config file: $config"
echo "Output dir: $output_dir"
echo "Risk table file: $risk_table_file"
echo "Barcode labels: $barcode_labels"

