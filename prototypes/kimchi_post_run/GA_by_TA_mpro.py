import pandas as pd
import os
import sys

if __name__ == "__main__":
    report_df = pd.read_csv(sys.argv[1], sep = "\t", error_bad_lines = False)
    
    ID_df = report_df["Taxonomy"].apply(lambda x: x.split(";")[-1]) #str.split(";", expand = True)
    #ID_df["reads"] = report_df["Reads"]
    #print(ID_df)
    #print(ID_df.columns())
    #6 = the partial taxa
    
    group_df = pd.DataFrame(ID_df)#[6].to_frame()
    group_df["reads"] = report_df["Reads"]
    group_df.columns = (["taxa", "reads"])
    #print(group_df)
    #group_df.to_csv("dummy.csv")
    group_df = group_df.groupby(["taxa"], as_index = False).sum()
    #sys.exit("breakier")
    sample_file = open(sys.argv[2], "r")
    sample_list = []
    for line in sample_file:
        sample_line = "s_"+ line
        sample_list.append(sample_line.strip("\n"))
    
    frames_list = []
    for item in sample_list:
        temp_df = group_df[group_df["taxa"].str.contains(item)]
        frames_list.append(temp_df)
    final_df = pd.concat(frames_list)    
        
    final_df.to_csv(sys.argv[3], index = False)
    #print("----------------------------")
    #print(group_df)
    #print(group_df[group_df["taxa"].str.contains("gasicomitatum")])
    