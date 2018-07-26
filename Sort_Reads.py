#!/usr/bin/env python
# this sorts the reads (they may come jumbled, or sparse), 
# IN: raw reads
# OUT: sorted reads (how? that's in the code)
import sys
import os
import os.path
import shutil
import subprocess
import multiprocessing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time

input_file = sys.argv[1]
output_file = sys.argv[2]

seqs = {}
start_import_time = time.clock()
with open(input_file, "r") as infile:
    seq_found = False
    seq_id = ""
    skip_plus = True
    line_count = 0
    for line in infile:
        if line.startswith("@") and not seq_found:
            seq_id = line[1:].split(" ")[0].strip("\n")
            seqs[seq_id] = ["",""]
            seq_found = True
        elif seq_found and seqs[seq_id][0] == "":
            seqs[seq_id][0] = line.strip("\n")
        elif line.startswith("+") and skip_plus:
            skip_plus = False
        elif seq_found and not skip_plus and seqs[seq_id][0] != "":
            seqs[seq_id][1] = line.strip("\n")
            seq_found = False
            skip_plus = True
        line_count += 1
        #if(line_count > 10):
        #    continue
            #sys.exit()
        #else:
        #    print(seq_id)
            
end_import_time = time.clock()
print("read file time:", end_import_time - start_import_time, "s")

def write_fastq(fastq_dict, mode):
    with open(output_file, mode) as outfile:
        for seq_id in fastq_dict:
            outfile.write("@" + seq_id + "\n")
            outfile.write(fastq_dict[seq_id][0] + "\n")
            outfile.write("+" + "\n")
            outfile.write(fastq_dict[seq_id][1] + "\n")

file_made = False
start_sort_time = time.clock()
sorted_seqs = sorted(seqs)
end_sort_time = time.clock()
print("sort time:", end_sort_time - start_sort_time, "s")
seq_records = {}

start_write_time = time.clock()
for seq_id in sorted_seqs:
    if seq_id.endswith(":1") or seq_id.endswith(":2") or seq_id.endswith(":3") or seq_id.endswith(":4"):
        seq_records[seq_id[:-2]] = seqs[seq_id]
    else:
        seq_records[seq_id] = seqs[seq_id]
    if len(seq_records) == 100000:
        if not file_made:
            file_made = True
            write_fastq(seq_records, "w")
        else:
            write_fastq(seq_records, "a")
        seq_records = {}
if len(seq_records) > 0:
    if not file_made:
        file_made = True
        write_fastq(seq_records, "w")
    else:
        write_fastq(seq_records, "a")
        
end_write_time = time.clock()
print("write time:", end_write_time - start_write_time, "s")
print("total time:", end_write_time - start_import_time, "s")