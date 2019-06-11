# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 18:40:29 2019

@author: Billy
"""

import os 
import sys
import pandas as pd

if __name__ == "__main__":
    mpro_taxa_df = pd.read_csv(sys.argv[1], sep = "\t").groupby(["Taxonomy"]).sum()
    mpro_taxa_df["Taxonomy"] = mpro_taxa_df.index
    mpro_taxa_df.reset_index(inplace = True, drop = True)
    
    print(mpro_taxa_df)