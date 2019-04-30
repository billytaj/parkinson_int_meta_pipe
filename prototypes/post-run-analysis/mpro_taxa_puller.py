import sys
import os
import pandas as pd


def make_genus_label (taxa):
    taxa_tree = taxa.split(";")
    if(len(taxa_tree) < 7):
        return "g_unknown"
    else:
        if(taxa_tree[6] == "g_"):
            return "g_unknown"
        else:
            return taxa_tree[6]



if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    
    selected_df = gene_df[~gene_df["Taxonomy"].str.contains("unclassified")]
    selected_df["RPK"] = 1
    selected_df["RPK"] = selected_df["RPK"].mask(selected_df["RPK"] > 0, selected_df["Reads"] / (selected_df["Length"] / 1000))
    print(selected_df)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found_old.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found.csv")
    selected_df = selected_df.groupby(["Taxonomy"], as_index = False).sum()
    
    
    
    selected_df["genus"] = selected_df["Taxonomy"].apply(lambda x: make_genus_label(x))
    print(selected_df)
    selected_df = selected_df.groupby(["genus"], as_index = False).sum()
    
    
    
    selected_df.to_csv(new_taxa_name, index = False)