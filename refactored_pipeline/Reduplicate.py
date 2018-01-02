#!/usr/bin/env python
#This file is just supposed to repopulate the mRNA with the right number of duplicates
#but in order to do so, it must find out what's mRNA.  
#There's better ways to do this, but we're told that it's such a minor thing, we should leave it be.
#It also sorta ruins the modularity aspect, since we now have to include more than the original 3 files, but.... yeah...
#We need a clever solution for this.
import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock

start_all = clock()
#Not actually called in this code
#Barrnap = "/home/j/jparkins/mobolaji/Tools/Barrnap/bin/barrnap"
#Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch"
#Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"

reference_file = sys.argv[1]
reference_sequences = SeqIO.to_dict(SeqIO.parse(reference_file, "fastq"))
dedeplicated_file = sys.argv[2]
dedeplicated_sequences = SeqIO.index(dedeplicated_file, "fastq")
cluster_file = sys.argv[3]
cluster_map = {}
reduplicated_file = sys.argv[4]
reduplicated_ids = set()
reduplicated_seqs = []

with open(cluster_file, "r") as clustr_read:
    rep = ""
    seq_id = ""
    for line in clustr_read:
        if line.startswith(">"):
            continue
        elif line.startswith("0"):
            rep = line[line.find(">") + 1:line.find("...")]
            seq_id = rep
            cluster_map[rep] = [seq_id]
        elif len(line) > 5:
            seq_id = line[line.find(">") + 1:line.find("...")]
            cluster_map[rep].append(seq_id)

for sequence in dedeplicated_sequences:
    if len(cluster_map[sequence]) > 1:
        for seq_id in cluster_map[sequence]:
            reduplicated_ids.add(seq_id)
    else:
        reduplicated_ids.add(sequence)

reduplicated_seqs = [reference_sequences[seq_id] for seq_id in sorted(reduplicated_ids)]

with open(reduplicated_file, "w") as out:
    SeqIO.write(reduplicated_seqs, out, "fastq")
    
end_all = clock()
print "Reduplicate"
print "================================"
print "total runtime:", end_all - start_all, "s"   