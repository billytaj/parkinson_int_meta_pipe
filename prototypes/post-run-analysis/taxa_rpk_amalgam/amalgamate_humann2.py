#this amalgamates the humann2 data
import sys
import os
import pandas as pd


if __name__ == "__main__":
    list_of_files = os.listdir(os.getcwd())
    df_count = 0
    base_df = None
    for item in list_of_files:
        if("humann2_taxa_found" in item):
            print("WORKING ON:", item)
            if (df_count == 0):
                df_count += 1
                base_df = pd.read_csv(item, index_col = 0)
                
                #print(base_df)
                #print("==================================================================================")
            else:
                print("WORKING ON:", item)
                new_df = pd.read_csv(item, index_col = 0)
                #print(new_df)
                #same_key_df = new_df(new_df["taxa"].isin(base_df["taxa"]))
                #diff_key_df = new_df(~new_df["taxa"].isin(base_df["taxa"]))
                base_df = base_df.add(new_df, fill_value = 0)
                df_count += 1
                #print(base_df)
                #print("==================================================================================")
    
    #print("DF count:", df_count)
    base_df["RPK"] = base_df["RPK"].mask(base_df["RPK"] > 0, base_df["RPK"]/df_count)
    total_rpk = base_df["RPK"].sum()
    base_df["percentage_rpk"] = 1
    base_df["percentage_rpk"] = base_df["percentage_rpk"].mask(base_df["percentage_rpk"] > 0, 100* base_df["RPK"] / total_rpk)
    base_df.sort_values("percentage_rpk", ascending = False, inplace = True)
    base_df.to_csv("humann2_avg_taxa_coverage.csv")