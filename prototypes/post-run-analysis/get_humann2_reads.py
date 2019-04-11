#This will read both excels, and get the read counts.

import sys
import os
import pandas as pd

def get_length_df(read_lengths_file):
    read_lengths_df = pd.read_csv(read_lengths_file, error_bad_lines = False)
    
    ID_df = read_lengths_df["GeneID"].str.split("|", expand = True)
    #print(ID_df)
    length_df = ID_df[[6, 7, 8]]
    print(length_df)
    length_df.columns = (["taxa", "uniref90", "uniref50"])
    length_df["length"] = read_lengths_df["Length"]
    return length_df

def get_rpk_df(gene_map_file):
    gene_map_df = pd.read_csv(gene_map_file, sep = "\t", error_bad_lines = False)
    gene_map_df.columns = ["GeneID", "RPK"]
    ID_df = gene_map_df["GeneID"].str.split("|", expand = True)
    ID_df["RPK"] = gene_map_df["RPK"]
    ID_df.dropna(axis = 0, inplace = True)
    ID_df.columns = (["uniref90", "taxa", "RPK"])
    return ID_df
    
if __name__ == "__main__":
    gene_map_file = sys.argv[1] #the ones with rpk -> gene familes
    read_lengths_file = sys.argv[2] #the bowtie file
    
    length_df = get_length_df(read_lengths_file)
    rpk_df = get_rpk_df(gene_map_file)
    
    rpk_df.to_csv("rpk_df.csv", mode = "w", index = False)
    length_df.to_csv("length_df.csv", mode = "w", index = False)
    
    
    #result_df = length_df[length_df["uniref90"].isin(rpk_df["uniref90"])]
    #result_df["RPK"] = rpk_df["RPK"]
    
    result_df = rpk_df
    result_df["length"] = length_df["length"]
    result_df.to_csv("humann2_final.csv")
    print(result_df)
    #result = pd.merge(length_df, rpk_df, how = "right", on = "uniref90")
    #result.drop("uniref50", 1, inplace = True)
    #result.drop("taxa", 1, inplace = True)
    
    #result.drop_duplicates("uniref90", inplace = True)
    #result.to_csv("merged.csv", mode = "w", index = False)
    #print(result)
    
   
    
    
    