import sys
import os
import pandas as pd


def make_genus_label (taxa):
    taxa_tree = taxa.split(";")
    print(taxa_tree)
    if(taxa == "unclassified"):
        return "unclassified"
    
    elif(taxa == "root"):
        return "unclassified"



if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    
    selected_df = gene_df#gene_df[~gene_df["Taxonomy"].str.contains("unclassified")]
    selected_df["RPK"] = 1
    selected_df["RPK"] = selected_df["RPK"].mask(selected_df["RPK"] > 0, selected_df["Reads"] / (selected_df["Length"] / 1000))
    selected_df.sort_values("RPK", inplace = True, ascending = False)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_0.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    selected_df = selected_df.groupby(["Taxonomy"], as_index = False).sum()
    selected_df.sort_values("RPK", inplace = True, ascending = False)
    print(selected_df)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_all.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    
    