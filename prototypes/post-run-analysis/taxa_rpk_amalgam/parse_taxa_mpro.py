# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:40:29 2019

@author: Billy
"""

import os 
import sys
import pandas as pd
from datetime import datetime as dt 
import multiprocessing as mp


def turn_to_list(item):
    final_list = []
    name_list = item.split(";")
    for line in name_list:
        final_item = line.split("_", 1)[-1]
        if not(final_item == ""):
            final_list.append(final_item)
    return final_list

def fetch_name(item):
    return item[-1]

def strip_prefix(item):
    if item:
        final_item = item.split("_", 1)[-1]
        return final_item
    else:
        return None
#def convert_to_taxid(item, names_df, nodes_df):
#    for entry in item:

def convert_name_to_taxid(name, names_df):
    #print(dt.today(), "working on:", name)
    if name is None:
        #print("Empty, return None")
        return None
    else:
            
        selected_df = names_df[names_df["name"].str.contains(name)].iloc[0]
        if not(selected_df.empty):
            taxa = selected_df["taxa"]
            #print("selected DF for key:", name, type(taxa), taxa)
            #print(selected_df)
            return taxa
        else:
            #print("EMPTY selected DF for key:", name)
            return 0
        
def convert_taxid_to_name(taxa_list, names_df):
    taxa_id = taxa_list[-1]
    selected_df = names_df[names_df["taxa"] == taxa_id].iloc[0]
    return selected_df["name"]
    
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
        #taxa_list.append(taxa_id)
        #taxa_dict[taxa_id] = taxa_category
        #count = 0
        #print("taxa id:", taxa_id)
        #while (~selected_nodes.empty):
        while(taxa_id != 1):    
            #print("working:", count)
            #count += 1
            selected_nodes = nodes_df.loc[nodes_df["taxa"] == taxa_id]
            
            
            taxa_category = selected_nodes["level"].iloc[0]
            #if(taxa_id == 131567):
            #    taxa_category = "Life on Earth"
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
    
        
def clean_list(item):
    cleaned_list = [x for x in item if str(x) != 'nan']
    cleaned_list = [x for x in cleaned_list if str(x) != 'None']
    cleaned_list = [x for x in cleaned_list if x != ""]
    if not cleaned_list:
        return None
    return cleaned_list        

def get_name(item):
    print("working on:", item)
    return item[-1]

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
        
        common_ancestor = 1
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
    

if __name__ == "__main__":
    
    manager = mp.Manager()
    return_dict = manager.dict()
    
    #sample_category_dict = manager.dict()
    sample_tree_return_dict = manager.dict()
    
    mpro_file = sys.argv[1]
    names_file = sys.argv[2]
    nodes_file = sys.argv[3]
    asf_list_file = sys.argv[4]
    output_location = sys.argv[5]
    
    
    mpro_taxa_df = pd.read_csv(mpro_file, sep = "\t").groupby(["Taxonomy"]).sum()
    mpro_taxa_df["Taxonomy"] = mpro_taxa_df.index
    mpro_taxa_df.reset_index(inplace = True, drop = True)
    #mpro_taxa_df["Taxonomy"] = mpro_taxa_df["Taxonomy"].apply(lambda x: turn_to_list(x))
    #mpro_taxa_df["name"] = mpro_taxa_df["Taxonomy"].apply(lambda x: fetch_name(x))
    
    taxa_df = mpro_taxa_df["Taxonomy"].str.split(";", expand = True)
    taxa_df = taxa_df.applymap(lambda x: strip_prefix(x))
    #print(taxa_df)
    #print(mpro_taxa_df)
    
    
    print(dt.today(), "started importing nodes and names")
    names_df = pd.read_csv(names_file, sep = "\t", header = None)
    nodes_df = pd.read_csv(nodes_file, sep = "\t", header = None)
    
    names_df = names_df.iloc[:, 0:8:2]
    names_df.columns = ["taxa", "name", "unique", "class"]
    nodes_df = nodes_df.iloc[:, 0:5:2]
    nodes_df.columns = ["taxa", "parent", "level"]
    
    print(dt.today(), "finished nodes and names import")
    
    print(dt.today(), "starting on name->number conversion")
    #taxa_df = taxa_df.applymap(lambda x: convert_name_to_taxid(x, names_df))
    print(dt.today(), "finished name->number conversion")
    
    taxa_df.columns = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxa_df.drop(["domain"], axis = 1, inplace = True)
    taxa_df["list_form"] = taxa_df.values.tolist()
    taxa_df["list_form"]= taxa_df["list_form"].apply(lambda x: clean_list(x))
    taxa_df = taxa_df.dropna(how = 'all')
    taxa_df.to_csv("tester.csv")
    print(taxa_df)
    taxa_df["name"] = taxa_df["list_form"].apply(lambda x: get_name(x))
    taxa_df = taxa_df[["name", "list_form"]]
    taxa_df["reads"] = mpro_taxa_df["Reads"]
    print(taxa_df)
    
    reads_df = taxa_df.groupby("name").sum()
    print("READS DF")
    print(reads_df)
    
    sample_list = taxa_df["name"].tolist()
    
    #sample_tree_dict, sample_category_dict = make_ref_dict(names_df, nodes_df, sample_list)
    

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
    
    
    
    
    #find common ancestry needs a list of numbers.  We have a list of names.  
    #we need to translate the list of names into a list of numbers.  
    
    print(dt.today(),"dealing with final results")
    final_dict = dict(return_dict)
    #taxa_df.index = taxa_df["name"]
    
    final_df = pd.Series(final_dict).to_frame("taxa_id")
    print(final_df)
    final_df["name"] = final_df.index
    final_df["reads"] = reads_df["reads"]
    final_df["genes"] = 1
    
    #print("READS DF")
    #print(reads_df)
    
    
    final_df["taxa_id"] = final_df["taxa_id"].apply(lambda x: ref_category_dict[x])
    inal_df = final_df.groupby("taxa_id", as_index = False).sum()
    
    final_df.to_csv(output_location + ".csv", index = False)
    print(final_df)
    
    #print(mpro_taxa_df)
    #print(taxa_df)