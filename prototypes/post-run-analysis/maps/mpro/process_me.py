#this script groups everything into genus-and-down, and looks for the ASF genus

import sys
import os
import pandas as pd

def rename_taxa(taxa):
    new_taxa = taxa
    #if(not taxa.endswith(";")):
    #    new_taxa = taxa + ";"
    if(taxa.endswith("\t")):
        new_taxa = taxa.rstrip()
    if("uncultured" in taxa):
        print(new_taxa)
        new_taxa = "unclassified"
    if("unidentified" in taxa):
        new_taxa = "unclassified"
    if(taxa == "root"):
        new_taxa = "unclassified"
    return new_taxa



if __name__ == "__main__":
    taxa_table_file = sys.argv[1]
    taxa_df = pd.read_csv(taxa_table_file, sep = "\t")
    
    taxa_df["Taxonomy"] = taxa_df["Taxonomy"].apply(lambda x: rename_taxa(x))
    taxa_split_df = taxa_df["Taxonomy"].str.split(";", expand = True)
    
    taxa_split_df.columns = ["root", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxa_split_df.fillna(value = "0", inplace = True)
    taxa_split_df["Reads"] = taxa_df["Reads"]
    #print(taxa_split_df)
    
    species_df = taxa_split_df[(taxa_split_df["species"] != "0")]#.groupby(["species"], as_index = False).sum()
    
    taxa_split_df = taxa_split_df[(taxa_split_df["species"] == "0")]
    genus_df = taxa_split_df[(taxa_split_df["genus"] != "0")]#.groupby(["genus"], as_index = False).sum()
    species_df.to_csv("demo_species.csv", index = False)
    print(species_df)