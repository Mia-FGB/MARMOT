import json
import pandas as pd

# Read the JSON file
with open("download.json", "r") as file:
    json_data = json.load(file)

# Convert JSON data to a pandas DataFrame
df_json = pd.DataFrame(json_data)

# Extract species from the JSON DataFrame
species_json = set(df_json['species'])


# Read the download.txt file
df_txt = pd.read_csv("download.txt", sep='\t', header=None, names=['species', 'taxid', 'assembly_accession', 'md5_url', 'fasta_url'])
# Extract species from the download.txt DataFrame
species_txt = set(df_txt['species'])

# Count the number of unique species in the JSON
species_count = df_json['species'].value_counts()
unique_species_count = df_json['species'].nunique()

# Find species that are in JSON but not in download.txt
species_only_in_json = species_json - species_txt

# Find species that are in download.txt but not in JSON
species_only_in_txt = species_txt - species_json

# # Print the count
# print("Number of species:", len(species_count))
# print("Number of unique species:", unique_species_count)

# Print the results
print("Species only in JSON:")
for species in species_only_in_json:
    print(species)

print("\nSpecies only in download.txt:")
for species in species_only_in_txt:
    print(species)