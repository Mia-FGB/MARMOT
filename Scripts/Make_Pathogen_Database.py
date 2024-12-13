import pandas as pd 
import json
import csv
import argparse
from ete3 import NCBITaxa
import requests

#add command line arguments ===
parser = argparse.ArgumentParser(description="Process PHIbase and Risk Register input files.")
parser.add_argument("-p", "--phibase", required=True, help="Path to the PHIbase input CSV file")
parser.add_argument("-r", "--risk_register", required=True, help="Path to the Risk Register input CSV file")

# Parse the command line arguments
args = parser.parse_args() 


# DEFRA Risk Register database ===
# https://planthealthportal.defra.gov.uk/pests-and-diseases/uk-plant-health-risk-register/downloadEntireRiskRegister.cfm 
#Download the most up to date version each time
#read in risk register
print("Reading in", args.risk_register)
risk_register = pd.read_csv(args.risk_register)
# Define the list of items to remove
remove = ["Insect", "Mite", "Nematode", "Plant"] 
# Filter rows where 'Type of pest' is not in the remove list
risk_register = risk_register[~risk_register['Type of pest'].isin(remove)] 
# Remove single quotes from 'Pest Name' column if present
risk_register['Pest Name'] = risk_register['Pest Name'].str.replace("'", "")


# PHIbase https://github.com/PHI-base/data/tree/master/releases ===
#Download the most up to date version each time 
#read in phibase 
print("Reading in", args.phibase)
phibase = pd.read_csv(args.phibase)

# Check for the correct host column (handle potential typo)
if "Host_description" in phibase.columns:
    host_column = "Host_description"
elif "Host_descripton" in phibase.columns:  # Note typo here
    host_column = "Host_descripton"
else:
    raise ValueError("Neither 'Host_description' nor 'Host_descripton' column found in PHIbase data.")

# Define the list of categories to keep
keep = ["monocots", "eudicots", "flowering plants",
        "basidiomycetes", "seed plants", "eukaryotes"]

# Filter the dataframe based on the host column values
phibase = phibase[phibase[host_column].isin(keep)]


#Merge the two datasets into one, only keeping species name ===
# Rename the columns to a common name ('species_name')
risk_register = risk_register[['Pest Name']].rename(columns={'Pest Name': 'species_name'})
phibase = phibase[['Pathogen_species']].rename(columns={'Pathogen_species': 'species_name'})

# Function to trim species names to only have the first three strings
def trim_species_name(name):
    return ' '.join(name.split()[:3])
# Apply the function to trim species names in both DataFrames
risk_register['species_name'] = risk_register['species_name'].apply(trim_species_name)
phibase['species_name'] = phibase['species_name'].apply(trim_species_name)

# Combine species names into a set to ensure uniqueness
unique_species = set(phibase['species_name']).union(set(risk_register['species_name']))
# Convert the set back to a dataframe
unique_species_df = pd.DataFrame({"species_name": list(unique_species)})
#note - could add extra species here if I wanted


# Initialize NCBI Taxa object ===
# this is instead of accessionTaxa.sql in R script
ncbi = NCBITaxa()
# Uncomment the following line only if you need to refresh the database (Downloads the taxonomy database)
#ncbi.update_taxonomy_database()


# Function to get TaxID for a species ===
def get_taxid(species_name):
    try:
        # Use NCBITaxa to fetch TaxID
        taxid = ncbi.get_name_translator([species_name])
        print("Getting TaxaID for", species_name)
        return taxid[species_name][0] if species_name in taxid else None
    except Exception as e:
        print(f"Error fetching TaxID for {species_name}: {e}")
        return None

# Add a TaxID column to the dataframe
unique_species_df['taxid'] = unique_species_df['species_name'].apply(get_taxid)
# Drop rows with missing TaxID
unique_species_df = unique_species_df.dropna(subset=['taxid'])

#Function to get Lineage from TaxID ===
# This is in the R script but I can't see why it is necessary as I have species names and taxids
# commenting out for now and if I don't need it will delete
# def get_lineage(taxid):
#     try:
#         # Fetch the lineage for the given TaxID
#         lineage = ncbi.get_lineage(taxid)
#         # Fetch the names corresponding to the lineage
#         names = ncbi.get_taxid_translator(lineage)
#         # Return the lineage names
#         return [names[taxid] for taxid in lineage]
#     except Exception as e:
#         print(f"Error fetching lineage for TaxID {taxid}: {e}")
#         return None

# # Add a lineage column to the dataframe
# unique_species_df['lineage'] = unique_species_df['taxid'].apply(get_lineage)

# Download Refseq & Genbank tables ===
# Function to download a file from a URL
def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check that the request was successful
    with open(filename, 'wb') as file:
        file.write(response.content)

# Download RefSeq and GenBank tables - Only need to do this occassionally to keep up to date
# print("Downloading RefSeq dataframe")
# download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt", "assembly_summary_refseq.txt")

# print("Downloading GenBank dataframe")
# download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", "assembly_summary_genbank.txt")

# Read in RefSeq dataframe
print("Reading in RefSeq dataframe")
refseq = pd.read_csv("assembly_summary_refseq.txt", sep='\t')
print(refseq.head())

# Example print statements to verify the input
# print(unique_species_df)
