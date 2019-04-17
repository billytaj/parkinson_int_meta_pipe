import os
import sys
import pandas as pd
import numpy as np

def add_ec(x):
    if(x == str(0)):
        return x
    else:   
        return "ec:" + x

if __name__ == "__main__":
    ec_map = sys.argv[1]
    humann2_ec_file = sys.argv[2]
    
    
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
    print(humann2_ec)
    #need to drop the no-ECs first
    
    
    metapro_matched = ec_df[ec_df["ec"].isin(humann2_ec["ec"])]
    print("metapro matched")
    print(metapro_matched)
    
    
    metapro_not_matched = ec_df[~ec_df["ec"].isin(humann2_ec["ec"])]
    print("metapro not matched")
    print(metapro_not_matched)
    
    humann2_not_matched = humann2_ec[~humann2_ec["ec"].isin(ec_df["ec"])]
    print("humann2 not matched")
    print(humann2_not_matched)