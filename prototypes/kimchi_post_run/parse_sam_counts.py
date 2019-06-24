import sys
import os
import pandas as pd
#use on the humann2 file
if __name__ == "__main__":
    counts_file = open(sys.argv[1], "r")
    counts_df = pd.read_csv(counts_file, sep = "\t", header = None)
    counts_df.columns = ["geneID", "reads"]
    print(counts_df.shape)
    id_df = counts_df["geneID"].str.rsplit("|", n = 4, expand = True)
    id_df["reads"] = counts_df["reads"]
    id_df.columns = ["the_rest","taxa", "uniref90", "uniref50",  "length", "reads"]
    id_df.drop(["the_rest"], axis = 1, inplace = True)
    id_df["genes"] = 1
    print(id_df)
    print(id_df.shape)
    #print(id_df)
    
    id_df = id_df.groupby(["taxa"], as_index = False).sum()
    print(id_df)
    counts_file.close()
    