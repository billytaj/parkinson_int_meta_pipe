#!/usr/bin/env python
#This file is just supposed to repopulate the mRNA with the right number of duplicates
#but in order to do so, it must find out what's mRNA.  
#There's better ways to do this, but we're told that it's such a minor thing, we should leave it be.
#It also sorta ruins the modularity aspect, since we now have to include more than the original 3 files, but.... yeah...
#We need a clever solution for this.
#Also, this needs to check for file existence.  If there's missing files in the arg, don't run
#also this thing is totally broken
import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock
import pandas as pd
start_all = clock()

def repopulate_double(ref_1_filename, ref_2_filename, mRNA_1_filename, mRNA_2_filename, cluster_1_filename, cluster_2_filename):
    #find the common ones between the 2 segments.  slot in whatever isn't common into orphans again.
    #not done
    ref_1_file = pd.read_csv(ref_1_filename, header = None, names = [None], sep = '\n', skip_blank_lines = False)
    ref_1_df = pd.DataFrame(ref_1_file.values.reshape(int(len(ref_1_file)/4), 4))
    ref_1_df.columns = ["ID", "seq", "junk", "quality"]

    ref_1_file = pd.read_csv(ref_2_filename, header = None, names = [None], sep = '\n', skip_blank_lines = False)
    ref_1_df = pd.DataFrame(ref_1_file.values.reshape(int(len(ref_1_file)/4), 4))
    ref_1_df.columns = ["ID", "seq", "junk", "quality"]
    
    
    
def repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename):
    ref_file = pd.read_csv(ref_filename, header = None, names = [None], sep = '\n', skip_blank_lines = False)
    ref_df = pd.DataFrame(ref_file.values.reshape(int(len(ref_file)/4), 4))
    ref_df.columns = ["ID", "seq", "junk", "quality"]
    #reference_sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fastq"))

    #mRNA_seqs = SeqIO.index(mRNA_file, "fastq")
    mRNA_file = pd.read_csv(mRNA_filename, header=None, names=[None], sep = '\n', skip_blank_lines = False)
    mRNA_df = pd.DataFrame(mRNA_file.values.reshape(int(len(mRNA_file)/4), 4))
    mRNA_df.columns = ["ID", "seq", "junk", "quality"]
    cluster_file = cluster_filename
    cluster_map = {}
    full_mRNA_file = output_filename

    #print("ref file:", reference_file)
    #print("unique file:", mRNA_file)
    #print("cluster file:", cluster_file)
    #print("full output:", full_mRNA_file)

    reduplicated_ids = set()
    reduplicated_seqs = []

    #There's some things that the duplicate remover deems as a duplicate, but comes as a different ID.  
    #scanning the cluster like this means we keep it all.  We must scan the cluster file
    with open(cluster_file, "r") as clustr_read:
        rep = ""
        seq_id = ""
        for line in clustr_read:
            if line.startswith(">"):
                continue
            elif line.startswith("0"):
                rep = "@" + line[line.find(">") + 1:line.find("...")]
                seq_id = rep
                cluster_map[rep] = [seq_id]
            elif len(line) > 5:
                seq_id = "@" + line[line.find(">") + 1:line.find("...")]
                cluster_map[rep].append(seq_id)
                
    
    #puts the cluster and mRNA IDs together in a single list
    for sequence in mRNA_df["ID"]:
        #print("looking at:", sequence)
        if(sequence in cluster_map):
            if len(cluster_map[sequence]) > 1:
                for seq_id in cluster_map[sequence]:
                    reduplicated_ids.add(seq_id)
                    #print("we're adding from the cluster map:", seq_id)
            else:
                reduplicated_ids.add(sequence)
                #print("adding from mRNA_df:", sequence)
        else:
            reduplicated_ids.add(sequence)
            #print("adding from mRNA_df:", sequence)

    #exports the full mRNA by fetching from ref 
    ref_df[ref_df.ID.isin(sorted(reduplicated_ids))].to_csv(full_mRNA_file, sep = '\n', mode = "w+", header = False, index = False)

    print("=================================")
    #print(ref_df[ref_df.ID.isin(sorted(reduplicated_ids))])

        
    end_all = clock()
    print ("Reduplicate")
    print ("================================")
    #print ("total runtime:", end_all - start_all, "s") 


if __name__ == "__main__":
    
    mRNA2_filename = None
    ref2_filename = None
    
    if(len(sys.argv) == 5):
        ref_filename = sys.argv[1]
        mRNA_filename = sys.argv[2]
        cluster_filename = sys.argv[3]
        output_filename = sys.argv[4]
        repopulate_single(ref_filename, mRNA_filename, cluster_filename, output_filename)
    
    elif(len(sys.argv) == 8):
        ref_filename = sys.argv[1]
        mRNA_filename = sys.argv[2]
        cluster_filename = sys.argv[3]
        output_filename = sys.argv[4]
        ref2_filename = sys.argv[5]
        mRNA2_filename = sys.argv[6]
        cluster2_filename = sys.argv[7]
    
    
  