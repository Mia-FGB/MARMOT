# 🐹 MARMOT

Metagenomic Alignment and Reporting for Monitoring Of Threats
MARMOT is a modular pipeline for detecting potential pathogens from Nanopore sequencing data. It performs quality filtering, alignment, and summary reporting — making it easy to go from raw reads to actionable microbial insights.

## Table of Contents
- [Usage](#usage)
- [Updating the Reference Database](#updating-the-reference-database)
- [Output Structure](#output-structure)
- [Scripts Overview](#scripts-overview)
- [Post-Processing with Taxonkit](#post-processing-with-taxonkit)

## Updating the Reference Database

The emergent pathogen database should be updated regularly. The latest update was on **08/05/24**.
The reference database can be updated using the following scripts, usually on local machine with interent access.
Conda environment - `pathogen_database`
**`build_reference_database.sh`**
These scripts can be found in https://github.com/Mia-FGB/MARMOT_ref_db/tree/main

## Usage

To use the pipeline, ensure the `config_template.sh` script is configured correctly with the following parameters:


## Config file
| Parameter             | Explanation                                                                                         | Example |
|-----------------------|-----------------------------------------------------------------------------------------------------|---------|
| `sample`              | Name of the sample being processed. Used for labeling and organization.                            | `CF_2023` |
| `location`            | Path to raw / rebasecalled reads. Should be the level above `fastq`.                               | `/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/raw/RL_24hCubTests_21092023/RL_24hCubTests_21092023/20230921_1707_X4_FAW72641_f261fc8c` |
| `filter_length`       | Reads shorter than this length will not be included in the analysis.                               | `300`   |
| `reference_database`  | Path to the reference database generated by `build_reference_database.sh`.                         | `/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/pathogen_database_080524.fa` |
| `scratch_dir`         | A directory for temporary large files to be stored. Can be the same as `output_dir`.               | `/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/pipeline_test` |
| `output_dir`          | Directory for output files.                                                                         | `/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/pipeline_test` |
| `barcode_list`        | List of barcodes to process. Can be a sequence or specific list.                                   | `$(seq -w 01 88)` or `("11" "02" "05" "08")` |
| `concatenated`        | Set to `"yes"` if the reads are already concatenated into one file, else `"no"`.                   | `yes`   |
| `contig_stats`        | Set to `"yes"` to calculate contig stats for the barcode, else `"no"`.                             | `yes`   |
| `genome_lengths_file` | File with TaxaIDs and genome lengths, generated by `build_reference_database.sh`.                  | `/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/080524_genome_lengths.tsv` |

---

### Submitting the Job Array

After configuring the file, submit your batch using the wrapper script:

```bash
./submit_with_config.sh config_template.sh
```

This wrapper will:
- Determine the number of barcodes from your config
- Create and submit a SLURM job array with the correct range (e.g. `--array=01-88`)
- Organize log files under `logs/<job_name>/`

## Output Structure

The pipeline generates the following output:
- **A directory for each barcode** containing result files and output/error messages for the individual jobs.
- **A log directory** storing output/error messages for each submission batch job.
- **Above the barcode directories:**
  - ` lcaparse_summary.txt`: Contains the Barcode, Read Count and % of reads assigned to each Taxa
  - Other summary statistics files.
  
## Scripts Overview

### `submit_batch_jobs.sh`
Runs `single_barcode_process.sh` for each barcode provided.

#### **Inputs:**
- Barcode numbers
- Raw read location
- Pre-filter length
- Reference database
- config

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
- `XX_num_fail.txt` (number of failed reads, only possible if there is a fastq_fail directory)
In Scratch_Dir
- `XX_barcode_300bp.fastq` (Reads filtered >300bp)

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


