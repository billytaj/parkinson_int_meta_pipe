import sys
import os
import pandas as pd


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
    
    
    
    taxa_df = selected_df["Taxonomy"].apply(lambda x: x.split(";"))
    print(taxa_df)
    
    
    selected_df.to_csv(new_taxa_name, index = False)