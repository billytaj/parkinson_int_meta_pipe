import os
import sys

import pandas as pd

#this script pulls out the ECs, matches against each metabolic pathway, and gives us a breakdown of how much each pathway is covered

#read in RPKM, and thepathway map


def turn_to_list(x):
    new_list = x.split("|")
    return new_list

if __name__ == "__main__":
    ec_path_file = sys.argv[1]
    rpkm_file = sys.argv[2]

    #get our ec path map searchable
    ec_path_df = pd.read_csv(ec_path_file)
    ec_path_df["Ecs"] = ec_path_df["Ecs"].apply(lambda x: turn_to_list(x))
    print(ec_path_df)
     
    num_ecs = ec_path_df["num_Ecs"]
    
    #get our annotated ECs
    
    rpkm_df = pd.read_csv(rpkm_file, sep = "\t")
    rpkm_df = rpkm_df[["GeneID", "EC#"]]
    rpkm_df = rpkm_df.drop_duplicates("EC#")
    rpkm_df = rpkm_df[rpkm_df["EC#"] != "0.0.0.0"]
    
    rpkm_ec = rpkm_df["EC#"].tolist()
    
    print(rpkm_ec)
    unique_ec = len(rpkm_ec)
    
    frames_list = []
    for item in rpkm_ec:
        selected_df = ec_path_df[ec_path_df["Ecs"].apply(lambda x: item in x)]
        selected_df["tally"] = 1
        frames_list.append(selected_df)    

    result = pd.concat(frames_list)
    
    print(result)
    
    print("number of unique ecs:", unique_ec)
    print("number of individual frames grabbed:", len(frames_list))
    count = 0
    for item in frames_list:
        print(count, item.shape)
        count += 1
    
    print("LOOKING FOR:", rpkm_ec[0])
    print(frames_list[0])
    
    result = result.groupby("Pathway", as_index = False).sum()
    result["num_Ecs"] = result["num_Ecs"].mask(result["tally"] > 1, result["num_Ecs"] / result["tally"])
    result["coverage"] = 1
    result["coverage"] = result["coverage"].mask(result["tally"] > 0, result["tally"] *100 / result["num_Ecs"])
    print(result)
    
    