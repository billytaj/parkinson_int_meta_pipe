import pandas as pd
import os
import sys

if __name__ == "__main__":
    rpkm_table_file = sys.argv[1]
    rpkm_df = pd.read_csv(rpkm_table_file, sep = "\t", error_bad_lines = False)
    gene_id_df = rpkm_df["GeneID"].str.split("|", 0, expand = True)
    gene_id_df.columns = ["gi", "num", "ref", "nc", "c", "num2", "taxa", "uniref90", "uniref50"]
    unique_taxa = gene_id_df["taxa"].unique()
    unique_taxa = list(filter(None, unique_taxa))
    #print(unique_taxa)
    
    
    unique_taxa_file = open("unique_taxa.txt", "w")
    unique_taxa_file.write("unique taxa\n")
    for item in unique_taxa:
        unique_taxa_file.write(item)
        unique_taxa_file.write("\n")
    unique_taxa_file.close()
    