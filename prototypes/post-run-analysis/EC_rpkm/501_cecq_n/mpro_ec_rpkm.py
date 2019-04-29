#this script takes finds the EC's rpkms by linking genes to ECs

import sys
import os
import pandas as pd

if __name__ == "__main__":
    gene_ec_file = sys.argv[1]
    gene_info_file = sys.argv[2]
    
    
    gene_info_df = pd.read_csv(gene_info_file, sep = '\t')
    #print(gene_info_df)
    
    
    gene_ec_df = pd.read_csv(gene_ec_file, header = None, sep = '\t')
    gene_ec_df.columns = ["GeneID", "EC"]
    #print(gene_ec_df)
    
    
    ec_only_df = gene_info_df[gene_info_df["GeneID"].isin(gene_ec_df["GeneID"])]
    ec_only_df["EC"] = gene_ec_df["EC"]
    print(ec_only_df)
    
    print("orig:", gene_info_df.shape)