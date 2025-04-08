# 🐹 MARMOT

Metagenomic Alignment and Reporting for Monitoring Of Threats
MARMOT is a modular pipeline for detecting potential pathogens from Nanopore sequencing data. It performs quality filtering, alignment, and summary reporting — making it easy to go from raw reads to actionable microbial insights.

## Table of Contents
- [Usage](#usage)
- [Updating the Reference Database](#updating-the-reference-database)
- [Output Structure](#output-structure)
- [Scripts Overview](#scripts-overview)
- [Post-Processing with Taxonkit](#post-processing-with-taxonkit)

## Usage

To use the pipeline, ensure the `config_template.sh` script is configured correctly with:
- **Raw Read Location**
 - Raw reads / rebasecalled on the HPC, level above fastq folders
- **Pre-filter length**
 - Reads shorter than this length will be removed before alignment (e.g. 300bp)
- **Reference database path**
 - Location of the reference database .fa file
- **Scratch Output Directory**
 -  Directory where concatenated reads can outpu to
- **Output Directory**
 - Where the output of the pipline will be written
 - It is recommended to run the script from here so the submission log is written to the same place
- **Barcode List**
 - Barcodes to analyse, can be as a sequence or list of specific barcodes 
 - e.g. barcode_list=("01" "02" "05 "08) or barcode_list=($(seq -w 01 88))


Then, run the script from your desired output directory:
```sh
/path/to/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/submit_batch_jobs.sh config_template.sh
```
Can use array=1-num_barcodes, to ensure nly the correct number of jobs are submitted

## Updating the Reference Database

The emergent pathogen database should be updated regularly. The latest update was on **08/05/24**.
The reference database can be updated using the following scripts, usually on local machine with interent access.
Conda environment - `pathogen_database`
**`Make_Pathogen_Database.py`** 
- Takes CSV files from DEFRA & PHIbase as an input
- Outputs a JSON file with the download links for the genomes and MD5 checksum files
**`download.py`**
 - Uses the Output JSON file to download the reference genomes and create the reference database


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


