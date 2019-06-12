# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:40:29 2019

@author: Billy
"""

import os 
import sys
import pandas as pd
from datetime import datetime as dt 

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
    selected_df = names_df[names_df["name"].str.contains(name)]
    return selected_df["taxa"]
    
        

if __name__ == "__main__":
    
    mpro_file = sys.argv[1]
    names_file = sys.argv[2]
    nodes_file = sys.argv[3]
    
    
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
    
    taxa_df = taxa_df.applymap(lambda x: convert_name_to_taxid(x, names_df))
    
    print(taxa_df)