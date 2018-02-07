#!/usr/bin/env python
#This is supposed to figure out what's leftover as pair 1, pair 2, contigs, and orphans
#It was originally named Map read contigs
# code seems to:
# 1) sift through the SAM file, and divides it into 2 sets:  stuff that's been mapped, and stuff that's not mapped
import sys
import os
import pandas as pd



def filter_contig_comsumption():


if __name__ == "__main__":
    pair_1_path = sys.argv[1]       #fastq
    pair_2_path = sys.argv[2]       #fastq
    orphans_path = sys.argv[3]      #fastq 
    pair_sam_path = sys.argv[4]     #sam 
    orphans_sam_path = sys.argv[5]   #sam
    output_path = sys.argv[6]       #just a location
    
    pair_1_df = pd.read_csv(pair_1_path, header=None, names=[None], sep = '\n', skip_blank_lines = False)
    pair_2_df = pd.read_csv(pair_2_path, header=None, names=[None], sep = '\n', skip_blank_lines = False)
    orphans_df = pd.read_csv(orphans_path, header=None, names=[None], sep = '\n', skip_blank_lines = False)
    
    pair_1_df.columns = ["ID", "sequence", "junk", "quality"]
    pair_2_df.columns = ["ID", "sequence", "junk", "quality"]
    orphans_df.columns = ["ID", "sequence", "junk", "quality"]
    
    pair_sam_df = pd.read_csv(pair_sam_path, header=None, sep="\t")
    pair_sam_df.iloc[:, 1] = pair_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])
    
    orphans_sam_df = pd.read_csv(orphans_sam_path, header=None, sep="\t")
    orphans_sam_df.iloc[:, 1] = orphans_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])

    #select mapped and unmapped slices from pair
    mapped_pair_sam_df = "@" + pair_sam_df.loc[pair_sam_df.iloc[:, 1] == "1"].iloc[:, 0]
    mapped_pair_sam_df = mapped_pair_sam_df.drop_duplicates()
    mapped_pair_sam_df.columns = ["ID"]
    
    unmapped_pair_sam_df = "@" + pair_sam_df.loc[pair_sam_df.iloc[:, 1] == "0"].iloc[:, 0]
    unmapped_pair_sam_df = unmapped_pair_sam_df.drop_duplicates()
    unmapped_pair_sam_df.columns = ["ID"]
    
    unmapped_pair_sam_df = unmapped_pair_sam_df[~unmapped_pair_sam_df.isin(mapped_pair_sam_df)]
    
    #select mapped and unmapped slices from orphans
    mapped_orphans_sam_df = orphans_sam_df.loc[orphans_sam_df.iloc[:, 1] == "1"].iloc[:, 0]
    mapped_orphans_sam_df = mapped_orphans_sam_df.drop_duplicates()
    mapped_orphans_sam_df.columns = ["ID"]
    
    unmapped_orphans_sam_df = orphans_sam_df.loc[orphans_sam_df.iloc[:, 1] == "0"].iloc[:, 0]
    unmapped_orphans_sam_df = unmapped_orphans_sam_df.drop_duplicates()
    unmapped_orphans_sam_df.columns = ["ID"]
    
    unmapped_orphans_sam_df = unmapped_orphans_sam_df[~unmapped_orphans_sam_df.isin(mapped_orphans_sam_df)]
    
    #------------------------------
    #write it
    pair_1_df[pair_1_df.ID.isin(unmapped_pair_sam_df)].to_csv(output_path+"pair_1.fastq", sep='\n', mode = "w+", header=False, index=False)
    pair_2_df[pair_2_df.ID.isin(unmapped_pair_sam_df)].to_csv(output_path+"pair_2.fastq", sep='\n', mode = "w+", header=False, index=False)
    orphans_df[orphans_df.ID.isin(unmapped_orphans_sam_df)].to_csv(output_path+"orphans.fastq", sep='\n', mode = "w+". header=False, index=False)
