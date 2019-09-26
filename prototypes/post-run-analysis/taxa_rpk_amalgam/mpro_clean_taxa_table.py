#This one's for mpro.  Feed it the taxa report.

import os 
import sys
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
    taxa_file = sys.argv[1]
    taxa_df = pd.read_csv(taxa_file, sep = "\t")
    taxa_df["Taxonomy"] = taxa_df["Taxonomy"].apply(lambda x: rename_taxa(x))
    taxa_split_df = taxa_df["Taxonomy"].str.split(";", expand = True)
    #taxa_df = taxa_df.drop(labels = '8', axis = 1)
    
    taxa_split_df.columns = ["root", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxa_split_df.fillna(value = "0", inplace = True)
    print(taxa_split_df)
    
    
    taxa_split_df["Length"] = taxa_df["Length"]
    taxa_split_df["Reads"] = taxa_df["Reads"]
    taxa_split_df["RPK"] = 0
    taxa_split_df["RPK"] = taxa_split_df["RPK"].mask(taxa_split_df["RPK"] == 0, taxa_split_df["Reads"] / (taxa_split_df["Length"] / 1000))
    #species_df = taxa_split_df.groupby(["species"], as_index = False).sum()
    #print("Species read sum:", species_df["Reads"].sum())
    #print(species_df)
    
    #genus_df = taxa_split_df.groupby(["genus"], as_index = False).sum()
    #print(genus_df)
    #print("genus read sum:", genus_df["Reads"].sum())
    
    species_df = taxa_split_df[(taxa_split_df["species"] != "0")].groupby(["species"], as_index = False).sum()
    print("Species read sum:", species_df["Reads"].sum(), "Species RPK sum:", species_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["species"] == "0")]
    
    genus_df = taxa_split_df[(taxa_split_df["genus"] != "0")].groupby(["genus"], as_index = False).sum()
    print("genus read sum:", genus_df["Reads"].sum(), "genus RPK sum:", genus_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["genus"] == "0")]
    
    family_df = taxa_split_df[taxa_split_df["family"] != "0"].groupby(["family"], as_index = False).sum()
    print("family read sum:", family_df["Reads"].sum(), "family RPK sum:", family_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["family"] == "0")]
    
    order_df = taxa_split_df[taxa_split_df["order"] != "0"].groupby(["order"], as_index = False).sum()
    print("order read sum:", order_df["Reads"].sum(), "order RPK sum:", order_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["order"] == "0")]
    
    class_df = taxa_split_df[taxa_split_df["class"] != "0"].groupby(["class"], as_index = False).sum()
    print("class read sum:", class_df["Reads"].sum(), "class RPK sum:", class_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["class"] == "0")]
    
    phylum_df = taxa_split_df[taxa_split_df["phylum"] != "0"].groupby(["phylum"], as_index = False).sum()
    print("phylum read sum:", phylum_df["Reads"].sum(), "phylum RPK sum:", phylum_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["phylum"] == "0")]
    
    kingdom_df = taxa_split_df[taxa_split_df["kingdom"] != "0"].groupby(["kingdom"], as_index = False).sum()
    print("kingdom read sum:", kingdom_df["Reads"].sum(), "kingdom RPK sum:", kingdom_df["RPK"].sum())
    taxa_split_df = taxa_split_df[(taxa_split_df["kingdom"] == "0")]
    
    unclassified_df = taxa_split_df.groupby(["root"], as_index = False).sum()
    print(unclassified_df)
    #print(taxa_split_df)
    
    
    #print(taxa_split_df)
    #print(taxa_split_df[taxa_split_df["genus"] == "g_"])
    #print(taxa_split_df[taxa_split_df["family"] == "f_"])
    #print(taxa_split_df[taxa_split_df["order"] == "o_"])
    #print(taxa_split_df[taxa_split_df["class"] == "c_"])
    #print(taxa_split_df[taxa_split_df["phylum"] == "p_"])
    #print(taxa_split_df[taxa_split_df["kingdom"] == "k_"])
    
    
    #print(taxa_df.columns)
    