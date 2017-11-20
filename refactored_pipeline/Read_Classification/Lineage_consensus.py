#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

target_rank = sys.argv[1]
read_cutoff = sys.argv[2]
read1_file = sys.argv[3]
#read1_sequences = SeqIO.index(read1_file, "fastq")
read2_file = sys.argv[4]
#read2_sequences = SeqIO.index(read2_file, "fastq")
Output_folder = sys.argv[5]
consensus_classification_file = sys.argv[6]
consensus_classification_summary_file = sys.argv[7]


nodes_file = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
nodes = {}

with open(nodes_file, "r") as infile:
    for line in infile:
        cols = line.split("\t|\t")
        taxid = cols[0]
        parent = cols[1]
        rank = cols[2]
        nodes[taxid] = (parent, rank)
    else:
        nodes["0"] = ("0", "unclassified")

names_file = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
names = {}

with open(names_file, "r") as infile:
    for line in infile:
        cols = line.split("\t|\t")
        taxid = cols[0]
        name = cols[1]
        if "scientific name" in cols[3]:
            names[taxid] = name

Ranks = ["unclassified", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]

read_classifications = []

for x in range(len(sys.argv[8:])):
    classification_file = sys.argv[x + 8]
    classification = {}
    with open(c, "r") as infile:
        for line in infile:
            cols = line.split("\t")
            read = cols[1]
            taxid = cols[2].strip()
            if taxid in nodes:
                while nodes[taxid][1] != target_rank:
                    if taxid == "1":
                        break
                    if nodes[taxid][1] in Ranks:
                        if Ranks.index(nodes[taxid][1]) > Ranks.index(target_rank):
                            taxid = nodes[taxid][0]
                        else:
                            break
                    else:
                        taxid = nodes[taxid][0]
                else:
                    if nodes[taxid][1] == target_rank:
                        classification[read] = names[taxid]
                    else:
                        classification[read] = "Unclassified"
    read_classifications.append(classification)

consensus_classifications = {}

for classification in read_classifications:
    for read in classification:
        for class_check in read_classifications:
            if read in class_check:
                if classification[read] == class_check[read]:
                    continue
                else:
                    consensus_classifications[read] = "Unclassified"
                    break
        else:
            consensus_classifications[read] = classification[read]

taxonomic_bins = {}

for read in consensus_classifications:
    if consensus_classifications[read] in taxonomic_bins:
        taxonomic_bins[consensus_classifications[read]].append(read)
    else:
        taxonomic_bins[consensus_classifications[read]] = [read]

try:
    len(taxonomic_bins["Unclassified"])
except:
    taxonomic_bins["Unclassified"] = []

small_taxa = []

for taxa in taxonomic_bins:
    if taxa == "Unclassified":
        continue
    if len(taxonomic_bins[taxa]) < len(consensus_classifications) * float(read_cutoff):
        taxonomic_bins["Unclassified"].extend(taxonomic_bins[taxa])
        small_taxa.append(taxa)

for taxa in small_taxa:
    del taxonomic_bins[taxa]

print len(consensus_classifications)

'''
with open(consensus_classification_file, "w") as outfile:
    for taxa in taxonomic_bins:
        for read in taxonomic_bins[taxa]:
            outfile.write(read + "\t" + taxa + "\n")

with open(consensus_classification_summary_file, "w") as outfile:
    percentages = {}
    for taxa in taxonomic_bins:
        percentages[taxa] = str(float(len(taxonomic_bins[taxa])) / float(len(consensus_classifications)))
    for sorted_taxa in sorted(percentages.items(), key=lambda x: x[1], reverse=True):
        taxa = sorted_taxa[0]
        outfile.write(taxa + "\t" + str(len(taxonomic_bins[taxa])) + "\t" + str(float(len(taxonomic_bins[taxa])) / float(len(consensus_classifications)) * 100) + "\n")
''''''
try:
    os.mkdir(Output_folder)
except:
    shutil.rmtree(Output_folder)
    os.mkdir(Output_folder)

for taxa in taxonomic_bins:
    seqs1 = [read1_sequences[seq_id] for seq_id in taxonomic_bins[taxa]]
    seqs2 = [read2_sequences[seq_id] for seq_id in taxonomic_bins[taxa]]
    with open(os.path.join(Output_folder, taxa + "_" + os.path.basename(read1_file)), "w") as out:
        SeqIO.write(seqs1, out, "fastq")
    with open(os.path.join(Output_folder, taxa + "_" + os.path.basename(read2_file)), "w") as out:
        SeqIO.write(seqs2, out, "fastq")
'''