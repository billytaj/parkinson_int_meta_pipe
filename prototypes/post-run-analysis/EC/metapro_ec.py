import os
import sys
import pandas as pd
import numpy as np

def add_ec(x):
    if(x == str(0)):
        return "none"
    else:   
        return "ec:" + x

if __name__ == "__main__":
    ec_map = sys.argv[1]
    humann2_ec_file = sys.argv[2]
    output_header = sys.argv[3]
    if not os.path.exists(output_header):
        os.mkdir(output_header)
    
    output_header = output_header + "/" + output_header
    ec_df = pd.read_csv(ec_map)
    ec_df = ec_df.T
    ec_df = ec_df.fillna(999)
    ec_list = ec_df.values.flatten()
    ec_list = [x for x in ec_list if x != 999]
    ec_df = pd.DataFrame(ec_list)
    ec_df.columns = (["ec"])
    
    print(ec_df)
    
    humann2_ec = pd.read_csv(humann2_ec_file)
    humann2_ec["ec"] = humann2_ec["ec"].apply(lambda x: add_ec(x))
    humann2_ec = humann2_ec[humann2_ec["ec"] != "none"]
    print(humann2_ec)
    #need to drop the no-ECs first
    
    
    metapro_matched = ec_df[ec_df["ec"].isin(humann2_ec["ec"])]
    output_name = output_header + "_metapro_humann2_match.csv"
    metapro_matched.to_csv(output_name, mode = "w", index = False, header = None)
    print("metapro matched")
    print(metapro_matched)
    
    
    metapro_not_matched = ec_df[~ec_df["ec"].isin(humann2_ec["ec"])]
    output_name = output_header + "_metapro_only.csv"
    metapro_not_matched.to_csv(output_name, mode = "w", index = False, header = None)
    
    print("metapro not matched")
    print(metapro_not_matched)
    
    humann2_not_matched = humann2_ec[~humann2_ec["ec"].isin(ec_df["ec"])]
    
    output_name = output_header + "_humann2_only.csv"
    humann2_not_matched.to_csv(output_name, mode = "w", index = False, header = None)
    print("humann2 not matched")
    print(humann2_not_matched)
    
    summary_table_name = output_header + "_summary_table.txt"
    summary_table = open(summary_table_name, "w")
    matched_line = "matched: " + str(metapro_matched.shape[0]) + "\n"
    humann2_line = "humann2 only: " + str(humann2_not_matched.shape[0]) + "\n"
    metapro_line = "metapro only: " + str(metapro_not_matched.shape[0]) + "\n"
    
    summary_table.write(matched_line)
    summary_table.write(humann2_line)
    summary_table.write(metapro_line)
    summary_table.close()