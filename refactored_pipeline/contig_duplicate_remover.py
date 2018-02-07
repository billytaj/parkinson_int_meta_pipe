#!/usr/bin/env python
#This is supposed to figure out what's leftover as pair 1, pair 2, contigs, and orphans
#It was originally named Map read contigs
# code seems to:
# 1) sift through the SAM file, and divides it into 2 sets:  stuff that's been mapped, and stuff that's not mapped
import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock
import sys
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
    mapped_pair_sam_df = pair_sam_df.loc[pair_sam_df.iloc[:, 1] == "1"].iloc[:, :1]
    mapped_pair_sam_df.columns = ["ID", "flag"]
    unmapped_pair_sam_df = pair_sam_df.loc[pair_sam_df.iloc[:, 1] == "0"].iloc[:, :1]
    unmapped_pair_sam_df.columns = ["ID", "flag"]
    
    #select mapped and unmapped slices from orphans
    mapped_orphans_sam_df = orphans_sam_df.loc[orphans_sam_df.iloc[:, 1] == "1"].iloc[:, :1]
    mapped_orphans_sam_df.columns = ["ID", "flag"]
    unmapped_orphans_sam_df = orphans_sam_df.loc[orphans_sam_df.iloc[:, 1] == "0"].iloc[:, :1]
    unmapped_orphans_sam_df.columns = ["ID", "flag"]
    pair_1_df[]
    
    
start_all = clock()
paired_file1 = sys.argv[1]
paired1_seqs = SeqIO.index(paired_file1, "fastq")
paired_file2 = sys.argv[2]
paired2_seqs = SeqIO.index(paired_file2, "fastq")
unpaired_file = sys.argv[3]
unpaired_seqs = SeqIO.index(unpaired_file, "fastq")
paired_sam_file = sys.argv[4]
unpaired_sam_file = sys.argv[5]
output_file = sys.argv[6]


contig2read_map = {}
unmapped_paired_reads = set()
unmapped_paired1_seqs = []
unmapped_paired2_seqs = []
unmapped_unpaired_reads = set()
unmapped_unpaired_seqs = []


def sam_filter(sam):
    #resets the flag column to show whether it's been mapped (0), or unmapped(1)
    sam_df = pd.read_csv(sys.argv[1], header=None, sep="\t")
    sam_df.iloc[:, 1] = sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8]) #convert the flag into binary, then extract the "mapped" flag at position 4
    return sam_df
"""
    with open(sam, "r") as samfile:
        for line in samfile:
            if not line.startswith("@") and len(line) > 1:
                line_parts = line.split("\t")
                flag = bin(int(line_parts[1]))[2:].zfill(11)
                if flag[8] == "0":
                    if line_parts[2] in contig2read_map:
                        if line_parts[0] not in contig2read_map[line_parts[2]]:
                            contig2read_map[line_parts[2]].append(line_parts[0])
                    else:
                        contig2read_map[line_parts[2]] = [line_parts[0]]
                elif flag[8] == "1":
                    unmapped.add(line_parts[0])
"""
pair_sam_df = sam_filter(pair_sam_path, unmapped_paired_reads)
orphan_sam_df = sam_filter(orphan_sam_path, unmapped_unpaired_reads)

for read in unmapped_paired_reads:
    unmapped_paired1_seqs.append(paired1_seqs[read])
    unmapped_paired2_seqs.append(paired2_seqs[read])
with open(os.path.splitext(paired_file1)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_paired1_seqs, out, "fastq")
with open(os.path.splitext(paired_file2)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_paired2_seqs, out, "fastq")

for read in unmapped_unpaired_reads:
    unmapped_unpaired_seqs.append(unpaired_seqs[read])
with open(os.path.splitext(unpaired_file)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_unpaired_seqs, out, "fastq")

with open(output_file, "w") as outfile:
    for contig in contig2read_map:
        outfile.write(contig + "\t" + str(len(contig2read_map[contig])))
        for read in contig2read_map[contig]:
            outfile.write("\t" + read)
        else:
            outfile.write("\n")

end_all = clock()
print "Map read contigs"
print "==========================="
print "total runtime:", end_all - start_all, "s"          