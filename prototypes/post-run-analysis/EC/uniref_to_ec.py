#this script takes humann2's uniref and turns it to EC numbers

import sys
import os
import pandas as pd


def translate_uniref(ec_dict, key):
    try:
        return ec_dict[key]
    except:
        return 0
if __name__ == "__main__":
    ec_map_file = sys.argv[1]
    humann2_gene_file = sys.argv[2]
    out_path = sys.argv[3]
    
    ec_dict = dict()
    
    ec_file = open(ec_map_file, "r")
    count = 0
    for item in ec_file:
        count += 1
        cleaned_item = item.strip("\n")
        item_list = cleaned_item.split("\t")
        ec_number = item_list[0]
        item_list.pop(0)
        for obj in item_list:
            ec_dict[str(obj)] = ec_number
    #ec_dict["UNMAPPED"] = 0
        
    #for item in ec_dict:
    #    print(item, ":", ec_dict[item])
    
    gene_map = pd.read_csv(humann2_gene_file, sep = "\t", error_bad_lines = False)
    #print(gene_map)
    
    ID_df = gene_map["# Gene Family"].str.split("|", expand = True)[0].drop_duplicates()
    
    ID_df = pd.DataFrame(ID_df)
    ID_df.columns = (["uniref"])
    ID_df["ec"] = ID_df["uniref"].apply(lambda x: translate_uniref(ec_dict, x))
    out_name = out_path + "_uniref_to_ec.csv"
    #ID_df = ID_df[ID_df.ec != 0]
    ID_df.to_csv(out_name, mode = "w", index = False)
    print(ID_df)
    
    
    