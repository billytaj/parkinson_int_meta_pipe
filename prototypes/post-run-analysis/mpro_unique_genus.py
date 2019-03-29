import pandas as pd
import os
import sys

if __name__ == "__main__":
    report_file = sys.argv[1]
    
    report_df = pd.read_csv(report_file, sep = "\t", error_bad_lines = False)
    
    ID_df = report_df["GeneID"].str.split("|", expand = True)
    ID_df["reads"] = report_df["Reads"]
    print(ID_df)
    #6 = the partial taxa
    group_df = ID_df[6].to_frame()
    group_df["reads"] = report_df["Reads"]
    group_df.columns = (["taxa", "reads"])
    print(group_df)
    group_df = group_df.groupby(["taxa"]).sum()
    print("----------------------------")
    print(group_df)
    group_df.to_csv("unique_species_taxa.csv", sep = ",", mode = "w")
    group_df["taxa"] = group_df.index
    
    print(group_df)
    
    group_df.reset_index(inplace = True, drop = True)
    
    
    
    
    genus_df = group_df["taxa"].str.split(".", expand = True)[0].to_frame()
    genus_df["reads"] = group_df["reads"]
    genus_df.columns = ["genus", "reads"]
    genus_df = genus_df.groupby(["genus"]).sum()
    genus_df.to_csv("unique_genus_taxa.csv", sep = ",", mode = "w")
    
    print(genus_df)
    