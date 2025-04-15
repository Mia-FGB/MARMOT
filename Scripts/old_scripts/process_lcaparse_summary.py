#!/usr/bin/python
# Doing this in R now

import pandas as pd
import sys


def extract_species(data):
    # Filter for Taxon_Rank == species
    species_df = data[data["Taxon_Rank"] == "species"].copy()
    # Extract the last entry in the Taxon_Path
    species_df["Species"] = species_df["Taxon_Path"].str.split(",").str[-1]

    return species_df

def merge_risk(risk, species_df):
    # Ensure the Species column exists in both DataFrames
    if "Species" not in species_df.columns:
        raise KeyError("The 'Species' column is missing in the species DataFrame.")
    if "Species" not in risk.columns:
        raise KeyError("The 'Species' column is missing in the risk DataFrame.")
    
    # Merge the risk table with the species DataFrame - currently leaving in the species that aren't in the risk table
    # Can code them as "Not in risk table" if needed
    merged_df = pd.merge(species_df, risk, on="Species", how="left")
    
    return merged_df

def write_merged_file(merged_df, output_file):
    # Write the merged DataFrame to a new file
    merged_df.to_csv(output_file, sep="\t", index=False)
    print(f"Merged data written to {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_lcaparse_summary.py <lcaparse_summary.txt> <risk_table.csv>")
        sys.exit(1)

    summary_file = sys.argv[1]  # First argument: summary file
    risk_table = sys.argv[2]  # Second argument: risk table

    # Read the summary file
    summary = pd.read_csv(summary_file, sep="\t", header=0)
    # Read the risk table
    risk = pd.read_csv(risk_table, sep=",", header=0)

    # Extract species from the summary
    species_df = extract_species(summary)

    # Merge the risk table with the species DataFrame
    merged_df = merge_risk(risk, species_df)

    # Write the merged DataFrame to a new file
    write_merged_file(merged_df, "merged_lcaparse_summary.txt")

if __name__ == "__main__":
    main()
