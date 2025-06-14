#!/bin/bash
# Check for the correct number of arguments
if [ "$#" -ne 9 ]; then
    echo "Usage single_barcode_process.sh: $0 <barcode_number> <location> <filter_length> <reference_database> <scratch_dir> <output_dir> <concatenated> <contig_stats> <genome_lengths_file>"
    exit 1
fi

# Extract the arguments - these will be provided by the submit_batch_jobs.sh script
barcode_number="$1"
location="$2"
filter_length="$3"
reference_database="$4"
scratch_dir="$5"
output_dir="$6"
concatenated="$7"
contig_stats="$8"
genome_lengths_file="$9"

# Check if the specified location exists
if [ ! -d "$location" ]; then
    echo "Error: Location '$location' does not exist."
    exit 1
fi

# Calculate the barcode directory path (output directory/barcode<barcode_number>)
barcode_dir="$output_dir/barcode${barcode_number}"
# Create barcode &  log directory if they don't exist
mkdir -p "$barcode_dir"
mkdir -p "$barcode_dir/logs"


# Execute the prep_reads.sh script with the barcode and location
JOBID1=$(sbatch --mem 1G \
    -p ei-short \
    -o "barcode${barcode_number}/logs/${barcode_number}_prepreads.out" \
    --error "barcode${barcode_number}/logs/${barcode_number}_prepreads.err" \
    --job-name="${barcode_number}_prepreads" \
    --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/prep_reads.sh $barcode_number $location $filter_length $scratch_dir $output_dir $concatenated $contig_stats" | awk '{print $NF}')

# Check if the job submission was successful
if [ -z "$JOBID1" ]; then
    echo "Error: Failed to submit the prep_reads job."
    exit 1
fi

echo "Submitted prep_reads job ($JOBID1)"

while true; do
    JOB1_STATUS=$(sacct -j $JOBID1 | awk 'NR==3 {print $6}')
    #echo "Job status for $JOBID1 is $JOB1_STATUS" # Uncomment for debugging
    if [[ "$JOB1_STATUS" == *COMPLETED* ]]; then
        echo "prep_reads has completed."
        break
    elif [[ "$JOB1_STATUS" == *FAILED* ]]; then
        echo "Error: The prep_reads job ($JOBID1) failed."
        exit 1
    elif [[ "$JOB1_STATUS" == *RUNNING* || "$JOB1_STATUS" == *PENDING* ]]; then
        echo "prep_reads still running/pending..."
        sleep 120
    fi
done

echo "The prep_reads job ($JOBID1) completed successfully."


# Minimap Job after JOB1 is finished - specific parameters -N 100, -c and -x map-ont
JOBID2=$(sbatch \
    --dependency=afterok:$JOBID1 \
    --mem 50G \
    -p ei-short \
    -o "$barcode_dir/logs/${barcode_number}_minimap.out" \
    --error "$barcode_dir/logs/${barcode_number}_minimap.err" \
    --job-name="${barcode_number}_minimap2" \
    --wrap "minimap2 -N 100 -c -x map-ont ${reference_database} \"$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq\" \\
    > \"$barcode_dir/${barcode_number}_mapped.paf\"" | awk '{print $NF}')

echo "Submitted minimap job ($JOBID2) will run when prep reads ($JOBID1) finishes"
echo "Minimap run with the following flags -N 100 -c -x map-ont, on the > ${filter_length}bp filtered reads"
echo "Barcode ${barcode_number} is being mapped to ${reference_database}"

# Check if the job submission was successful
if [ -z "$JOBID2" ]; then
    echo "Error: Failed to submit the minimap job."
    exit 1
fi


while true; do
    JOB2_STATUS=$(sacct -j $JOBID2 | awk 'NR==3 {print $6}')
    #echo "Job status for $JOBID2 is $JOB2_STATUS" # Uncomment for debugging
    if [[ "$JOB2_STATUS" == *COMPLETED* ]]; then
        echo "minimap Job has completed."
        break
    elif [[ "$JOB2_STATUS" == *FAILED* ]]; then
        echo "Error: The minimap job ($JOBID2) failed. Killing overall job."
        exit 1
    elif [[ "$JOB2_STATUS" == *RUNNING* || "$JOB2_STATUS" == *PENDING* ]]; then
        echo "minimap Job still running/pending..."
        sleep 120
    fi
done

echo "The minimap job (${JOBID2}) completed successfully."

#Submit job for lca_parse script
JOBID3=$(sbatch \
    --dependency=afterok:$JOBID2 \
    --mem 75G \
    -p ei-short \
    -o "$barcode_dir/logs/${barcode_number}_lcaparse.out" \
    --error "$barcode_dir/logs/${barcode_number}_lcaparse.err" \
    --job-name="${barcode_number}_lcaparse" \
    --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/run_lcaparse.sh $barcode_number" | awk '{print $NF}')

# Check if the job submission was successful
if [ -z "$JOBID3" ]; then
    echo "Error: Failed to submit the lca_parse job."
    exit 1
fi

echo "Submitted lca_parse ($JOBID3) will run when minimap ($JOBID2) finishes"
echo "Barcode ${barcode_number} minimap paf file is being lca parsed"

while true; do
    JOB3_STATUS=$(sacct -j $JOBID3 | awk 'NR==3 {print $6}')
    #echo "Job status for $JOBID3 is $JOB3_STATUS" # Uncomment for debugging
    if [[ "$JOB3_STATUS" == *COMPLETED* ]]; then
        echo "lca_parse job has completed."
        break
    elif [[ "$JOB3_STATUS" == *FAILED* ]]; then
        echo "Error: The lca_parse job ($JOBID3) failed."
        exit 1
    elif [[ "$JOB3_STATUS" == *RUNNING* || "$JOB3_STATUS" == *PENDING* ]]; then
        echo "lca_parse Job still running/pending..."
        sleep 120
    fi
done

echo "The lca_parse job ($JOBID3) completed successfully."

# Genome coverage script
JOBID4=$(sbatch --dependency=afterok:$JOBID3 \
    -p ei-short \
    --mem 5G \
    -o "$barcode_dir/logs/${barcode_number}_genome_coverage.out" \
    -e "$barcode_dir/logs/${barcode_number}_genome_coverage.err" \
    --job-name="${barcode_number}_genome_coverage" \
    --wrap "python /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/Scripts/pathogen_genome_coverage_from_paf.py $barcode_number $genome_lengths_file $barcode_dir" | awk '{print $NF}')

# make sure the genome coverage job has completed and created the files before completing single_barcode_process.sh
while true; do
    JOB4_STATUS=$(sacct -j $JOBID4 | awk 'NR==3 {print $6}')
    #echo "Job status for $JOBID4 is $JOB4_STATUS" # Uncomment for debugging
    if [[ "$JOB4_STATUS" == *COMPLETED* ]]; then
        echo "genome coverage job has completed."
        break
    elif [[ "$JOB4_STATUS" == *FAILED* ]]; then
        echo "Error: The genome coverage job ($JOBID4) failed."
        exit 1
    elif [[ "$JOB4_STATUS" == *RUNNING* || "$JOB4_STATUS" == *PENDING* ]]; then
        echo "genome coverage job still running/pending..."
        sleep 120
    fi
done

echo "The genome coverage job ($JOBID4) completed successfully."

echo "Barcode ${barcode_number} has been processed"

