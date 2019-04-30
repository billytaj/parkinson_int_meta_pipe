import sys
import os
import pandas as pd


if __name__ == "__main__":
    gene_file = sys.argv[1]
    new_name = sys.argv[2]
    gene_df = pd.read_csv(gene_file, sep = "\t",  error_bad_lines = False)
    gene_df.columns = (['GeneID', 'RPK'])
    
    selected_df = gene_df[gene_df["GeneID"].str.contains("g__")]
    selected_df["taxa"] = selected_df["GeneID"].apply(lambda x:x.split("|")[1])
    new_taxa_name = os.path.join(new_name, new_name + "_humann2_taxa_found_old.csv")
    selected_df.to_csv(new_taxa_name, index = False)
    selected_df = selected_df.groupby(["taxa"], as_index = False).sum()
    print(selected_df)
    new_taxa_name = os.path.join(new_name, new_name + "_humann2_taxa_found.csv")
    selected_df.to_csv(new_taxa_name, index = False)