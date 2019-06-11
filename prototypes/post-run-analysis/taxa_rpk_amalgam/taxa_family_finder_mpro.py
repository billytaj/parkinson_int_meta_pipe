#this takes in taxa from 2 sides, and tries to find the family.
#then it's gotta find match up with the reference sets, and where.

import sys
import os
import pandas as pd
import multiprocessing as mp
from datetime import datetime as dt

def construct_tree_num(names_df, nodes_df, num):
    
    #takes in a taxa number, and then constructs the tree from that.
    #searches until it can't find anything (root).  everything ends at 1
    
    taxa_list = []
    #taxa_category_list = []
    taxa_dict = dict()
        
    selected_nodes = nodes_df.loc[nodes_df["taxa"] == num]
    if not(selected_nodes.empty):
        taxa_id = num
        taxa_category = selected_nodes["level"].iloc[0]
        #taxa_category_list.append(taxa_category)
        taxa_list.append(taxa_id)
        taxa_dict[taxa_id] = taxa_category
        #count = 0
        while(taxa_id != 1):    
            #print("working:", count)
            #count += 1
            selected_nodes = nodes_df.loc[nodes_df["taxa"] == taxa_id]
            
            taxa_category = selected_nodes["level"].iloc[0]
            if(taxa_id == 131567):
                taxa_category = "Life on Earth"
            if not (taxa_category == "no rank"):
            #taxa_category_list.append(taxa_category)
                taxa_dict[taxa_id] = taxa_category
                taxa_list.append(taxa_id)
                
            taxa_id = selected_nodes["parent"].iloc[0]
            
    taxa_list.reverse()
    return taxa_list, taxa_dict
    #return taxa_list, taxa_category_list


def construct_tree_name(names_df, nodes_df, name):
    
    selected_df = names_df[names_df["name"].str.contains(name)]
    taxa_list = []
    #taxa_category_list = []
    taxa_dict = dict()    
    
    selected_nodes = nodes_df.loc[nodes_df["taxa"].isin(selected_df["taxa"])]
    if not(selected_nodes.empty):
        #taxa_tree.append(selected_df["taxa"])
        taxa_id = selected_df["taxa"].iloc[0]
        taxa_category = selected_nodes["level"].iloc[0]
        #taxa_category_list.append(selected_nodes["level"].iloc[0])
        taxa_list.append(taxa_id)
        taxa_dict[taxa_id] = taxa_category
        #count = 0
        #print("taxa id:", taxa_id)
        #while (~selected_nodes.empty):
        while(taxa_id != 1):    
            #print("working:", count)
            #count += 1
            selected_nodes = nodes_df.loc[nodes_df["taxa"] == taxa_id]
            
            
            taxa_category = selected_nodes["level"].iloc[0]
            if(taxa_id == 131567):
                taxa_category = "Life on Earth"
            if not(taxa_category == "no rank"):
                taxa_dict[taxa_id] = taxa_category
                taxa_list.append(taxa_id)
            
            taxa_id = selected_nodes["parent"].iloc[0]
        #print("-=-=-=-=-=-=-=-=-=-=-=-=-=")
        #print("taxa list")
        #print(taxa_list)
        
        #print("-=-=-=-=-=-=-=-=-=-=-=")
        #print("taxa category")
        #print(taxa_category_list)
    taxa_list.reverse()
    return taxa_list, taxa_dict    
    #return taxa_list, taxa_category_list

#def check_common_ancestry(ref_dict, sample_list)
    

def call_me(names_df, nodes_df, name):
    #taxa_list, taxa_category_list = construct_tree_name(names_df, nodes_df, name)
    taxa_list, taxa_dict = construct_tree_name(names_df, nodes_df, name)
    print("================================")
    print(name)
    print(taxa_list)
    #print(taxa_category_list)
    print(taxa_dict)
    
    

def call_me_maybe(names_df, nodes_df, num):
    #taxa_list, taxa_category_list = construct_tree_num(names_df, nodes_df, num)
    taxa_list, taxa_dict = construct_tree_num(names_df, nodes_df, num)
    print("--------------------------------")
    print(num)
    print(taxa_list)
    #print(taxa_category_list)
    print(taxa_dict)





def make_ref_dict(names_df, nodes_df, asf_list):
    ref_tree_dict = dict()
    ref_category_dict = dict()    
    for item in asf_list:
        taxa_list, taxa_dict = construct_tree_name(names_df, nodes_df, item)
        ref_tree_dict[item] = taxa_list
        ref_category_dict = {**ref_category_dict, **taxa_dict}
    return ref_tree_dict, ref_category_dict


def make_sample_tree_dict(names_df, nodes_df, sample_name, tree_dict):
    #we don't actually need the categories from the sample, since they'll be matching the ones from the ASF anyway
    taxa_list, taxa_dict = construct_tree_name(names_df, nodes_df, sample_name)
    tree_dict[sample_name] = taxa_list
    


def find_common_ancestry(sample_name, ref_dict, sample_list, return_dict):
    #This will search through all of ASF, and figure out which matching is the highest, and return that.`
    prior_depth = 0
    for ref_key in ref_dict:
        
        ref_list = ref_dict[ref_key]
        ref_size = len(ref_list)
        sample_size = len(sample_list)
        range_of_loop = 0
        if(ref_size > sample_size):
            range_of_loop = sample_size
        else:
            range_of_loop = ref_size
        
        common_ancestor = -1
        search_depth = 0
        for i in range(0, range_of_loop):
            search_depth += 1
            if(ref_list[i] == sample_list[i]):
                common_ancestor = sample_list[i]
            else:
                break
        #print(sample_name, "last ancestor against", ref_key, "->", )
        
        if(search_depth > prior_depth):        
            return_dict[sample_name] = common_ancestor
            prior_depth = search_depth

def import_samsa_file(samsa_file):
    samsa_df = pd.read_csv(samsa_file, sep = "\t", header = None)
    samsa_df.columns = ["percent", "reads", "name"]
    samsa_list = samsa_df["name"].tolist()
    return samsa_list

def import_mpro_file(mpro_file):
    mpro_df = pd.read_csv(mpro_file, sep = "\t")
    

def export_from_dict(final_dict, ref_category_dict):
    if not final_dict:
        print("return dict is empty.  something's wrong")
    else:
        print("return dict not empty.  proceed")
        
        with open("samsa_out.txt", "w") as out:
            
            for key in final_dict:
                last_match_taxa = final_dict[key]
                category = "none"
                if last_match_taxa in ref_category_dict:
                    category = ref_category_dict[last_match_taxa]
                else:
                    category = "not found"
                out.write(key + " correct up to: " + category + "\n")

    #print(sample_name, "last common ancestor:", return_dict[sample_name])
if __name__ == "__main__":
    
    manager = mp.Manager()
    return_dict = manager.dict()
    
    #sample_category_dict = manager.dict()
    sample_tree_return_dict = manager.dict()
    
    names_file = sys.argv[1]
    nodes_file = sys.argv[2]
    asf_list_file = sys.argv[3]
    samsa_file = sys.argv[4]
    output_name = sys.argv[5]
    
    print(dt.today(), "started importing nodes and names")
    names_df = pd.read_csv(names_file, sep = "\t", header = None)
    nodes_df = pd.read_csv(nodes_file, sep = "\t", header = None)
    
    names_df = names_df.iloc[:, 0:8:2]
    names_df.columns = ["taxa", "name", "unique", "class"]
    nodes_df = nodes_df.iloc[:, 0:5:2]
    nodes_df.columns = ["taxa", "parent", "level"]
    #print(nodes_df)
    
    print(dt.today(), "finished importing nodes and names")
    #call_me(names_df, nodes_df, "sp. KRMCY2")
    
    #prep the references
    print(dt.today(), "started crafting ASF tree")
    asf_list = []
    with open(asf_list_file) as asf_file:
        for line in asf_file:
            asf_list.append(line.rstrip("\n"))
    
    ref_tree_dict, ref_category_dict = make_ref_dict(names_df, nodes_df, asf_list)
    
    for item in ref_tree_dict:
        print(item, ref_tree_dict[item])
    print(dt.today(), "finished crafting ASF tree")
        
        
    #call_me(names_df, nodes_df, "Parabacteroides goldsteinii")
    
    #print("=====================")
    #print("375288")
    #print(nodes_df[nodes_df["taxa"] == 375288])
    
    """
    #prep the samsa
    print(dt.today(), "started crafting sample tree dict")
    sample_list = import_samsa_file(samsa_file)
    
    print(dt.today(), "starting sample tree maker jobs")
    #sample_tree_dict, sample_category_dict = make_ref_dict(names_df, nodes_df, samsa_list)
    sample_tree_maker_job_store = []
    for item in sample_list:
        sample_tree_maker_job = mp.Process(
            target = make_sample_tree_dict,
            args = (names_df, nodes_df, item, sample_tree_return_dict)
            )
        sample_tree_maker_job_store.append(sample_tree_maker_job)
        sample_tree_maker_job.start()
    print(dt.today(), "sample tree jobs launched")
    for item in sample_tree_maker_job_store:
        item.join()
    sample_tree_maker_job_store[:] = []
    sample_tree_dict = dict(sample_tree_return_dict)
    
    print(dt.today(), "finished crafting sample tree dict")
    
    

    mp_store = []
    print(dt.today(), "started launching jobs")
    for sample_key in sample_tree_dict:
        sample_list = sample_tree_dict[sample_key]
        mp_job = mp.Process(
            target = find_common_ancestry, 
            args = (sample_key, ref_tree_dict, sample_list, return_dict)
            )
        mp_store.append(mp_job)
        mp_job.start()
    print(dt.today(), "finished launching jobs")    
    for item in mp_store:
        item.join()
    mp_store[:] = []
    print(dt.today(), "finished running jobs")
    
    """
    
    
    #find common ancestry needs a list of numbers.  We have a list of names.  
    #we need to translate the list of names into a list of numbers.  
    
    print(dt.today(),"dealing with final results")
    final_dict = dict(return_dict)
    
    samsa_df.index = samsa_df["name"]
    
    final_df = pd.Series(final_dict).to_frame("taxa_id")
    final_df["name"] = final_df.index
    final_df["reads"] = samsa_df["reads"]
    final_df["genes"] = 1
    
    print(final_df)
    final_df["taxa_id"] = final_df["taxa_id"].apply(lambda x: ref_category_dict[x])
    final_df = final_df.groupby("taxa_id", as_index = False).sum()
    
    final_df.to_csv(output_name + ".csv", index = False)
    print(final_df)
    
    #print("+++++++++++++++++++++++++++++++++")
    #print(samsa_df)
    
    #export_from_dict(final_dict, ref_category_dict)
    
