import sys
import os
import pandas as pd

if __name__ == "__main__":
    gene_names_file = sys.argv[1]
    output_name = sys.argv[2] #just give it some 
    gene_names_df = pd.read_csv(gene_names_file, sep = "\t")
    gene_names_df = gene_names_df.groupby("English name").sum()
    gene_names_df.drop(["read length"], axis = 1, inplace = True)
    gene_names_df.to_csv(output_name + "_summary.csv")
    print(gene_names_df)