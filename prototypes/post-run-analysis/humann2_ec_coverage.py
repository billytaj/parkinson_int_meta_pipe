#this pulls out all ECs from Humann2, the ec -> path map, and tells us what the coverage level is for each pathway

import pandas as pd
import sys
import os

def turn_to_list(x):
    new_list = x.split("|")
    return new_list

if __name__ == "__main__":
    ec_path_file = sys.argv[1]
    humann2_file = sys.argv[2]

    #get our ec path map searchable
    ec_path_df = pd.read_csv(ec_path_file)
    ec_path_df["Ecs"] = ec_path_df["Ecs"].apply(lambda x: turn_to_list(x))
    print(ec_path_df)
     
    num_ecs = ec_path_df["num_Ecs"]
    
    humann2_df = pd.read_csv(humann2_file)
    humann2_df = humann2_df[humann2_df["ec"] != "0"]
    humann2_df = humann2_df.drop_duplicates("ec")
    humann2_ec = humann2_df["ec"].tolist()
    
    frames_list = []
    for item in humann2_ec:
        selected_df = ec_path_df[ec_path_df["Ecs"].apply(lambda x: item in x)]
        selected_df["tally"] = 1
        frames_list.append(selected_df)    

    result = pd.concat(frames_list)
    
    result = result.groupby("Pathway", as_index = False).sum()
    result["num_Ecs"] = result["num_Ecs"].mask(result["tally"] > 1, result["num_Ecs"] / result["tally"])
    result["coverage"] = 1
    result["coverage"] = result["coverage"].mask(result["tally"] > 0, result["tally"] *100 / result["num_Ecs"])
    print(result)