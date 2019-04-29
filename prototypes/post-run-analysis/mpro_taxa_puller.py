import sys
import os
import pandas as pd


if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    
    selected_df = gene_df[~gene_df["Taxonomy"].str.contains("unclassified")]
    new_taxa_name = os.path.join(new_name, new_name + "_mpro_taxa_found.csv")
    selected_df.to_csv(new_taxa_name, index = False)