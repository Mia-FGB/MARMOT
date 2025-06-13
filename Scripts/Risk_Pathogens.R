# Packages 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readr)
library(fs)  # for file path safety
library(stringr)
library(RColorBrewer)
library(viridis)

cat("Commencing R script... \n")

# Set paths and read in the data -----
# Script starts the same as create_risk_plots.R
# Won't include graph plotting, will have more parameter options for risk data

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Expecting: Rscript create_risk_plots.R <output_dir> <risk_table> <barcode_labels>
if (length(args) < 3) {
  stop("Usage: Rscript create_risk_plots.R <output_dir> <risk_table> <barcode_labels>")
}

# output_dir is where the write file output is saved
output_dir <- args[1]
risk_table <- args[2]
barcode_labels <- args[3]

# Check if the input directory exists
if (!dir.exists(output_dir)) {
  stop("Output directory does not exist: ", output_dir)
}
# Check if the risk table file exists & read in
if (!file.exists(risk_table)) {
  stop("Risk table file does not exist: ", risk_table)
}
cat("Reading in risk table ", risk_table, "\n")
risk_table <- read_csv(risk_table,  show_col_types = FALSE)


# Check if the barcode labels file exists & read in
if (!file.exists(barcode_labels)) {
  stop("Barcode labels file does not exist: ", barcode_labels, )
}
cat("Reading in barcode labels ", barcode_labels, "\n")

barcode_labels <- 
  read_delim(barcode_labels, delim = "\t",  show_col_types = FALSE)

barcode_labels$Barcode <- as.character(barcode_labels$Barcode)
barcode_labels$Barcode_Label <- as.character(barcode_labels$Barcode_Label)

cat("Barcode label file read successfully:\n")
cat("  Number of rows: ", nrow(barcode_labels), "\n")
cat("  Unique Barcodes in labels:\n")
print(head(unique(barcode_labels$Barcode)))

cat("Set paths and read in the data, output_dir ", output_dir, "\n")
# Construct output directory path
graph_save_path <- path(output_dir, "Graphs")

# Check and create the graph output directory if it doesn't exist
if (!dir.exists(graph_save_path)) {
  dir.create(graph_save_path, recursive = TRUE)
}

# Function to read input files
read_input_file <- function(filename, output_dir, delim = "\t") {
  readr::read_delim(fs::path(output_dir, filename),
  delim = delim, show_col_types = FALSE)
}

# Read all required files, these were created in the write_files.sh script
lcaparse_perread <- read_input_file("lcaparse_perread.txt", output_dir)
genome_coverage <- read_input_file("genome_coverage_all.txt", output_dir)
read_numbers     <- read_input_file("read_numbers.tsv", output_dir)

# Convert Barcode columns to character for consistency
lcaparse_perread$Barcode <- as.character(lcaparse_perread$Barcode)
read_numbers$Barcode <- as.character(read_numbers$Barcode)
genome_coverage$Barcode <- as.character(genome_coverage$Barcode)



# Set plotting theme
custom_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    axis.line.y = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90)
  )


# Process data ------------------------------------------------------------

# Risk data ------
# Filter out viruses as their species names aren't correct
risk_table <- risk_table %>%
  filter(Type_of_pest != ("Virus or Viroid"))

# Look at the risk data - see how many species are present more than once
risk_counts <- risk_table %>%
  count(Species, sort = TRUE) %>%
  filter(n > 1)

# Candidatus Phytoplasma is in there 25 times, then 11 species present 2 - 6 times

# Collapsing to have one species per row, taking the max value of each
# If there are 5 strains of a Bacteria we are assuming worst case scenario
collapsed_risk_table <- risk_table %>%
  group_by(Species) %>%
  summarise(
    Regulation = paste(
      unique(na.omit(Regulation)), 
      collapse = ";"
    ),
    
    Likelihood_unmitigated = if (all(is.na(Likelihood_unmitigated))) NA_real_ 
      else max(Likelihood_unmitigated, na.rm = TRUE),
    Impact_unmitigated = if (all(is.na(Impact_unmitigated))) NA_real_ 
      else max(Impact_unmitigated, na.rm = TRUE),
    Risk_Rating_unmitigated = if (all(is.na(Risk_Rating_unmitigated))) NA_real_ 
      else max(Risk_Rating_unmitigated, na.rm = TRUE),
    
    Likelihood_mitigated = if (all(is.na(Likelihood_mitigated))) NA_real_ 
      else max(Likelihood_mitigated, na.rm = TRUE),
    Impact_mitigated = if (all(is.na(Impact_mitigated))) NA_real_ 
      else max(Impact_mitigated, na.rm = TRUE),
    Risk_Rating_mitigated = if (all(is.na(Risk_Rating_mitigated))) NA_real_ 
      else max(Risk_Rating_mitigated, na.rm = TRUE),
    
    Regulated = if_else(
      any(Regulated == "Yes", na.rm = TRUE), 
      "Yes", 
      "No"
    ),
    Natural_Spread = if_else(
      any(Natural_Spread == "Yes", na.rm = TRUE), 
      "Yes", 
      "No"
    ),
    
    UK = if (length(unique(na.omit(UK))) == 1) 
      unique(na.omit(UK)) 
      else "N/A",
    
    .groups = "drop"
  )


#Lcaparse data ----

# Make anything below species level species instead e.g. subspecies to species 
# This does not handle species group / complex
lcaparse_perread <- lcaparse_perread %>%
  mutate(
    Species = case_when(
      Taxon_Rank == "subspecies" ~ str_extract(Taxon_Name, "^[^ ]+\\s[^ ]+(-[^ ]+)?"),
      TRUE ~ Taxon_Name
    )
  )

# Summarise and filter the file to have one row per species per barcode
# Lose the taxonID in this step as it is complicated by the combined subspecies 
lcaparse <- lcaparse_perread %>%
  group_by(Barcode, Species) %>%
  summarise(
    Taxon_Name = first(Taxon_Name),
    Taxon_Rank = first(Taxon_Rank),
    Read_Count = n_distinct(Read_ID),
    Taxon_ID = first(Taxon_ID),
    Avg_Mean_Identity = mean(Mean_Identity, na.rm = TRUE),
    .groups = "drop"
  )

# merge read numbers & lcaparse on Barcode
read_numbers <- read_numbers %>% 
  rename(TotalReadCount = Read_Count)

lcaparse <- lcaparse %>% 
  left_join(read_numbers, by = "Barcode")

#Create normalised read count columns 
lcaparse <- lcaparse %>% 
  mutate(
    HP100k = (Read_Count / TotalReadCount) * 100000,
    Filtered_HP100k = (Read_Count / FilterReadCount) * 100000)


# Combining Read & DEFRA Risk Data    -----
# This will loose the unassigned reads so if I want them in the graph would have to add back in later
lca_risk <- lcaparse %>%
  left_join(collapsed_risk_table, by = "Species") %>% 
  filter(Read_Count >10) # Filter low level species that could be artefacts, may need to play with filter 

# Add a column that groups by risk factor - groupings based on DEFRA documentation
lca_risk  <- lca_risk  %>%
  mutate(
    Risk_Category = case_when(
      Risk_Rating_mitigated <= 14                     ~ "Blue",
      Risk_Rating_mitigated >= 15 & Risk_Rating_mitigated <= 29  ~ "Green",
      Risk_Rating_mitigated >= 30 & Risk_Rating_mitigated <= 44  ~ "Yellow",
      Risk_Rating_mitigated >= 45 & Risk_Rating_mitigated <= 59  ~ "Orange",
      Risk_Rating_mitigated >= 60                              ~ "Red",
      TRUE ~ "Unclassified"
    )
  )

# Improve UK column - unknown to be n/a
lca_risk$UK[is.na(lca_risk$UK) | lca_risk$UK == "Unknown"] <- "N/A"


lca_risk$Risk_Category <- factor(
  lca_risk$Risk_Category,
  levels = c("Red", "Orange", "Yellow", "Green", "Blue", "Unclassified")
)

#For debugging 
cat("Joining barcode labels to lca_risk...\n")
cat("  Unique Barcodes in lca_risk before join:\n")
print(head(unique(lca_risk$Barcode)))

# Add Barcode labels to the lca_risk data
lca_risk <- lca_risk %>%
  left_join(barcode_labels, by = "Barcode")

# Replace Barcode column with Label for plotting
lca_risk$Barcode_Label <- factor(
  lca_risk$Barcode_Label,
  levels = unique(barcode_labels$Barcode_Label)  # preserve input order
)

# For debugging
cat("After join:\n")
cat("  Total rows in lca_risk: ", nrow(lca_risk), "\n")
cat("  Rows with missing Barcode_Label: ", sum(is.na(lca_risk$Barcode_Label)), "\n")
cat("  Example Barcode_Label values:\n")
print(head(unique(lca_risk$Barcode_Label)))

missing_labels <- setdiff(unique(lca_risk$Barcode), barcode_labels$Barcode)
cat("Barcodes in lca_risk missing from barcode_labels:\n")
print(missing_labels)

#Summarise the data to only include species with a DEFRA defined Risk category
risk_only <- lca_risk %>% 
  filter(Risk_Category != ("Unclassified"))



# Pull out the species and maybe ReadIDs to consider in more depth -----
extract_red_risk_reads <- function(data_risk,
                                data_perread,
                                output_dir,
                                risk_category,
                                exclude_widespread = FALSE,
                                filename = "RedRisk_ReadIDs.tsv") {
# Base filtering (not on read count or identity)
filtered <- data_risk %>%
filter(
    Risk_Category == risk_category,
)

# Optional additional filter
if (exclude_widespread) {
filtered <- filtered %>% filter(UK != "Present (Widespread)")
}

# Get read IDs with metadata
barcode_taxa <- filtered %>%
distinct(Barcode, Taxon_Name)

read_ids <- data_perread %>%
semi_join(barcode_taxa, by = c("Barcode", "Taxon_Name")) %>%
select(Barcode, Read_ID, Taxon_Name) %>%
left_join(
    filtered %>% select(Barcode, Taxon_Name, Taxon_ID, Risk_Category, UK, Avg_Mean_Identity, Read_Count),
    by = c("Barcode", "Taxon_Name")
)

# Save file
save_path <- fs::path(output_dir, filename)
write_tsv(read_ids, save_path)

cat("Saved ", nrow(read_ids), "read IDs for", risk_category, "risk to:" , save_path, "\n")
}

# Red risk not considering presence in the UK
extract_red_risk_reads(
data_risk = risk_only,
data_perread = lcaparse_perread,
output_dir = output_dir,
risk_category = "Red",
exclude_widespread = FALSE,
filename = "RedRisk_ReadIDs_all.tsv"
)

# Red Risk, filtering out those already known to be widespread present in the UK
extract_red_risk_reads(
data_risk = risk_only,
data_perread = lcaparse_perread,
output_dir = output_dir,
risk_category = "Red",
exclude_widespread = TRUE,
filename = "RedRisk_ReadIDs_noWidespread.tsv"
)

# For Orange Risk species -----
extract_red_risk_reads(
data_risk = risk_only,
data_perread = lcaparse_perread,
output_dir = output_dir,
risk_category = "Orange",
exclude_widespread = FALSE,
filename = "OrangeRisk_ReadIDs_all.tsv"
)

# Red Risk, filtering out those already known to be widespread present in the UK
extract_red_risk_reads(
data_risk = risk_only,
data_perread = lcaparse_perread,
output_dir = output_dir,
risk_category = "Orange",
exclude_widespread = TRUE,
filename = "OrangeRisk_ReadIDs_noWidespread.tsv"
)

# Create a summary file of species of interest, for all risk_categroies -----
# No longer filtering on Read count or Avg_Mean_Identity
create_risk_summary <- function(data_risk,
                                genome_coverage,
                                output_dir,
                                risk_category = NULL,
                                filename) {

  # Apply filter conditionally
  if (!is.null(risk_category)) {
    filtered_data <- data_risk %>%
      filter(Risk_Category == risk_category)
  } else {
    filtered_data <- data_risk
  }
  
  summary_table <- filtered_data %>%
    select(
      Barcode, Barcode_Label, Taxon_Name,
      Taxon_Rank, Taxon_ID, Risk_Category,
      UK, Read_Count, Avg_Mean_Identity
    ) %>%
    distinct() %>%
    left_join(
      genome_coverage %>%
        select(Barcode, taxaID, coverage_percentage, num_reads),
      by = c("Barcode" = "Barcode", "Taxon_ID" = "taxaID")
    ) %>%
    rename(
      Genome_Coverage = coverage_percentage,
      Coverage_Reads = num_reads
    )
  
  # Save
  save_path <- fs::path(output_dir, filename)
  write_tsv(summary_table, save_path)
  cat("Saved", risk_category, "risk summary table to:", save_path, "\n")
}

# Create summaries for all risk categories
risk_categories <- c("Red", "Orange", "Yellow", "Green", "Blue")
for (cat in risk_categories) {
  create_risk_summary(
    data_risk = risk_only,
    genome_coverage = genome_coverage,
    output_dir = output_dir,
    risk_category = cat,
    filename = paste0(cat, "_Risk_Species_Summary.tsv")
  )
}

create_risk_summary(
  data_risk = risk_only,
  genome_coverage = genome_coverage,
  output_dir = output_dir,
  risk_category = NULL,
  filename = "All_Risk_Species_Summary.tsv"
)

# Create Taxon summary file  -----
create_taxon_presence_summary <- function(data_risk,
                                          output_dir,
                                          filename = "Taxon_Presence_Summary.tsv") {
  total_samples <- n_distinct(data_risk$Barcode)                                      
  summary_table <- data_risk %>%
    group_by(Taxon_Name) %>%
    summarise(
      Risk_Category = paste(unique(Risk_Category)),
      UK = paste(unique(UK)),
      Num_Samples_Detected = n_distinct(Barcode),
      Proportion_Samples_Detected = (Num_Samples_Detected / total_samples) * 100,
      Max_Read_Count = max(Read_Count, na.rm = TRUE),
      Max_Identity = max(Avg_Mean_Identity, na.rm = TRUE),
      Avg_Read_Count = mean(Read_Count, na.rm = TRUE),
      Avg_Identity = mean(Avg_Mean_Identity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Num_Samples_Detected), desc(Max_Read_Count))
  
  save_path <- fs::path(output_dir, filename)
  write_tsv(summary_table, save_path)
  cat("Saved per-taxon presence summary to:", save_path, "\n")
}

create_taxon_presence_summary(
  data_risk = risk_only, 
  output_dir = output_dir
)

# To finish 
cat("R script for pathogens complete \n")
