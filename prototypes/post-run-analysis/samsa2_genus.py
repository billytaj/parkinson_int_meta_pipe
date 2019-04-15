#this attempts to sum up the genus levels of samsa2's outputs.  

import os
import sys
import pandas as pd
import numpy as np

if __name__ == "__main__":
    samsa2_file = sys.argv[1]
    #output_file = sys.argv[2]
    samsa2_df = None
    if samsa2_file.endswith(".tsv"):
        samsa2_df = pd.read_csv(samsa2_file, error_bad_lines = False, sep = "\t", header = None)
    elif(samsa2_file.endswith(".csv")):
        samsa2_df = pd.read_csv(samsa2_file, error_bad_lines = False, sep = ",", header = None)
    samsa2_df.columns = ["percent", "reads", "taxa"]
    
    
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Clostridium")) | (samsa2_df["taxa"].str.contains("ASF356")) | (samsa2_df["taxa"].str.contains("ASF502"))]
    selected_df.append(pd.Series([np.nan]), ignore_index = True)
    print("Clostridium")
    print(selected_df)
    print("Clostridium SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_clostridium.csv", mode = "w", index = False, header = None)
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Eubacterium")) | (samsa2_df["taxa"].str.contains("ASF492"))]
    selected_df.append(pd.Series([np.nan]), ignore_index = True)
    print("Eubacterium")
    print(selected_df)
    print("Eubacterium SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_eubacterium.csv", mode = "w", index = False, header = None)
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Firmicutes")) | (samsa2_df["taxa"].str.contains("ASF500"))]
    selected_df.append(pd.Series([np.nan]), ignore_index = True)
    print("Firmicutes")
    print(selected_df)
    print("Firmicutes SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_firmicutes.csv", mode = "w", index = False, header = None)
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Lactobacillus")) | (samsa2_df["taxa"].str.contains("ASF360")) | (samsa2_df["taxa"].str.contains("ASF361"))]
    selected_df.append(pd.Series([np.nan]), ignore_index = True)
    print("Lactobacillus")
    print(selected_df)
    print("Lactobacillus SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_lactobacillus.csv", mode = "w", index = False, header = None)
    
    
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Mucispirillum")) | (samsa2_df["taxa"].str.contains("ASF457"))]
    selected_df.append(pd.Series([np.nan]), ignore_index = True)
    print("Mucispirillum")
    print(selected_df)
    print("Mucispirillum SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_mucispirillium.csv", mode = "w", index = False, header = None)
    
    
    
    selected_df = samsa2_df.loc[(samsa2_df["taxa"].str.contains("Parabacteroides")) | (samsa2_df["taxa"].str.contains("ASF519"))]
    print("Parabacteroides")
    print(selected_df)
    print("Parabacteroides SUM:", selected_df["reads"].sum())
    #selected_df.to_csv(output_file + "_parabacteroides.csv", mode = "w", index = False, header = None)