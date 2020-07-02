#This filters our genes that didn't quite make it after the conversion of contig-segments to reads made them disappear.
#it'll remove the genes that didn't make it.

#use the final_gene_map

import sys
import pandas as pd
            
if __name__ == "__main__":
    EC_file = sys.argv[1]       #protein_ECs.all
    gene_map_file = sys.argv[2] #final_gene_map.tsv
    export_ec_file = sys.argv[3]
    
    
    EC_df = pd.read_csv(EC_file, sep = "\t", header = None, names = ["gene", "EC_count", "EC"])
    gene_map_df = pd.read_csv(gene_map_file, sep = "\t", usecols = [0, 1, 2])
    
    
    print(gene_map_df)
    print(EC_df)
    
    EC_df = EC_df[EC_df["gene"].isin(gene_map_df["gene"])]
    
    EC_df.to_csv(export_ec_file, sep = "\t", header = None, index = None)