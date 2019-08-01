import pandas as pd
import os
import sys
from datetime import datetime as dt
import multiprocessing as mp

def get_tally(p_id, return_dict, samsa_list, samsa_df, names_df, nodes_df):
    tally_dict = dict()
    for item in samsa_list:
        tally_key = find_taxa_rank(item, names_df, nodes_df)
        if(tally_key in tally_dict):
            current_read_count = tally_dict[tally_key]
            fetch_read_count = samsa_df[samsa_df["name"] == item]["reads"].iloc[0]
            #print("=================")
            #print("old key:", tally_key)
            #print("old count:", current_read_count)
            #print("updated fetch read", fetch_read_count)
            tally_dict[tally_key] = current_read_count + fetch_read_count
        else:
            fetch_read_count = samsa_df[samsa_df["name"] == item]["reads"].iloc[0]
            #print("-==-=-=-=-=-=-=-=-=-=-=-=-=-=-")
            #print("new key:", tally_key)
            #print("new fetch read:", fetch_read_count)
            tally_dict[tally_key] = fetch_read_count
        #print(tally_dict)
    return_dict[p_id] = tally_dict

def find_taxa_rank(name, names_df, nodes_df):
    selected_df = names_df[names_df["name"].str.contains(name)]
    #print("selected df:", selected_df)
    if(selected_df.empty):
        if("sp." in name):
            print(dt.today(), name, ": name mangled by SAMSA2.  it's likely a species")
            return "species"
        else:
            print(dt.today(), name, ": name mangled by SAMSA2.  it's unrecognizable")
            return "no rank"
    
    else:
        selected_taxa = selected_df["taxa"].iloc[0]
        #print("Selected taxa:", selected_taxa)
        selected_node = nodes_df[nodes_df["taxa"] == selected_taxa]
        #print("=================================")
        #print(name)
        #print("selected node:", selected_node)
        selected_rank = selected_node["level"].iloc[0]
        #print(selected_taxa, selected_rank)
        return selected_rank

if __name__ == "__main__":
    #gotta sort and figure out which name here is which level of taxa.  Then we can group.
    
    
    names_file = sys.argv[1]
    nodes_file = sys.argv[2]
    samsa_file = sys.argv[3]
    output_file = sys.argv[4]
    names_df = pd.read_csv(names_file, sep = "\t", header = None)
    names_df = names_df.iloc[:, 0:8:2]
    names_df.columns = ["taxa", "name", "unique", "class"]
    
    nodes_df = pd.read_csv(nodes_file, sep = "\t", header = None)
    nodes_df = nodes_df.iloc[:, 0:5:2]
    nodes_df.columns = ["taxa", "parent", "level"]
    
    
    samsa_df = pd.read_csv(samsa_file, sep = "\t", header = None)
    samsa_df.columns = ["percent", "reads", "name"]
    samsa_list = samsa_df["name"].tolist()
    
    print(names_df)
    
    print(nodes_df)
    print(dt.today(), "started apply")
    #samsa_df["rank"] = samsa_df["name"].apply(lambda x: find_taxa_rank(x, names_df, nodes_df))
    print(dt.today(), "started groupby")
    #samsa_df = samsa_df.groupby("rank", as_index = False).sum()
    
    print(dt.today(), "launching jobs")
    manager = mp.Manager()
    return_dict = manager.dict()
    p_id_count = 0
    process_list = []
    chunksize = 100
    for i in range(0, len(samsa_list), chunksize):
    
        list_segment = samsa_list[i:i+chunksize]
        tally_maker_process = mp.Process(
            target = get_tally,
            args = (p_id_count, return_dict, list_segment, samsa_df, names_df, nodes_df)
        )
        tally_maker_process.start()
        process_list.append(tally_maker_process)
        print(dt.today(), "launched job:", p_id_count)
        p_id_count += 1
        
    for item in process_list:
        item.join()
    process_list[:] = []
    
    print(dt.today(), "all jobs launched")
    final_dict = {}
    
    for item in return_dict:
        tally_dict = return_dict[item]
        for keys in tally_dict:
            if(keys in final_dict):
                cur_count = final_dict[keys]
                new_count = cur_count + tally_dict[keys]
                final_dict[keys] = new_count
            else:
                final_dict[keys] = tally_dict[keys]
        
        
    
    print(final_dict)
    #print(samsa_df)
    #samsa_df.to_csv(output_file + ".csv")
    # for item in samsa_list:
        # print("======================================")
        # print("SEARCHING for:", item)
        # selected_df = names_df[names_df["name"].str.contains(item)]
        # selected_df = selected_df.iloc[0]
        # #print("selected df:", selected_df)
        
        # selected_taxa = selected_df["taxa"]
        # #print("Selected taxa:", selected_taxa)
        # selected_node = nodes_df[nodes_df["taxa"] == selected_taxa]
        # #print("selected node:", selected_node)
        # selected_rank = selected_node["level"].iloc[0]
        # #print("selected rank:", selected_rank)
        # #print(selected_taxa, selected_rank)
        
        # print("SELECTED RANK:", str(selected_rank))