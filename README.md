# Emergent Pathogen pipeline scripts 

This repository contains scripts for detecting emergent pathogens in nanopore sequencing data. The pipeline processes raw nanopore reads, applies quality filtering, maps reads to a reference database, and generates statistics for downstream analysis.
Linked to - ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts

## Table of Contents
- [Usage](#usage)
- [Updating the Reference Database](#updating-the-reference-database)
- [Output Structure](#output-structure)
- [Scripts Overview](#scripts-overview)
- [Post-Processing with Taxonkit](#post-processing-with-taxonkit)

## Usage

To use the pipeline, ensure the `submit_batch_jobs.sh` script is configured correctly with:
- **Barcode numbers**
- **Pre-filter length**
- **Reference database path**
- **Raw read path**

Then, run the script from your desired output directory:
```sh
/path/to/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/submit_batch_jobs.sh
```

## Updating the Reference Database

The emergent pathogen database should be updated regularly. The latest update was on **08/05/24**.

## Output Structure

The pipeline generates the following output:
- **A directory for each barcode** containing result files.
- **A log directory** storing output/error messages for each batch job.
- **Above the barcode directories:**
  - `all_taxaID_counts.tsv`: Contains statistics for each barcode.
  - Other summary statistics files.
  
## Scripts Overview

### `submit_batch_jobs.sh`
Runs `single_barcode_process.sh` for each barcode provided.

#### **Inputs:**
- Barcode numbers
- Raw read location
- Pre-filter length
- Reference database

#### **Outputs:**
- Unzipped fastq file of all pass reads
- Unzipped fastq file of all reads > filter length

---
### `single_barcode_process.sh`
Submits a sequence of sbatch jobs in the following order:
1. `prep_reads.sh`
2. `minimap2`
3. `paf_parse.py`

#### **Inputs:**
Provided via `submit_batch_jobs.sh`:
- Barcode number
- Location
- Filter length

#### **Outputs:**
- All outputs from `prep_reads.sh`, `minimap2`, and `paf_parse.py`
- `percent_retained_results.txt`: Percentage of reads retained after filtering
- `fail_number_all.txt`: Number of failed reads per barcode
- `number_ignored_reads_from_parse.txt`: Reads ignored due to multiple mapping or barcode issues

---
### `prep_reads.sh`
Unzips and concatenates pass reads, generates contig statistics, and applies length filtering.

#### **Inputs:**
- Barcode
- Raw read location
- Filter length (e.g., 300)

#### **Outputs:**
Located in `barcodeXX` directory:
- `XX_contig_stats.txt`
- `XX_barcode.fastq` (all pass reads)
- `XX_num_fail.txt` (number of failed reads)
- `XX_barcode_300bp.fastq`

---
### `paf_parse.py`
Processes mapped `.paf` files, filtering out low-quality mappings (MQ < 5) and handling multiple mappings intelligently.

#### **Inputs:**
- `XX_mapped.paf`

#### **Outputs:**
- `XX_taxaID_count.txt` (tab-delimited TaxaID & counts)
- `XX_ignored_reads_query_ID.txt` (query, MQ values, multiple Taxa IDs)
- `XX_number_ignored_reads.txt` (number of ignored reads)

## Post-Processing with Taxonkit

Once all jobs have finished, post-processing is required to obtain lineage information:

1. Download `all_taxaID_counts.tsv` to your local machine.
2. Activate the `taxonkit` Conda environment:
   ```sh
   conda activate taxonkit
   ```
3. Use `taxonkit` to retrieve lineage information, outputting a CSV file.
4. This processed file is useful for **Phyloseq analysis**.

---
### **License & Citation**
Please cite this repository if you use it in your research. License details can be found in the repository.

For any issues or contributions, please submit a GitHub issue or pull request.


