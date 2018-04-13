#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
from time import clock as clock


input_file = sys.argv[1]
input_seqs = SeqIO.to_dict(SeqIO.parse(input_file, "fastq"))
BLAT_tab_file = sys.argv[2]
output_file = sys.argv[3]
output_seqs = []
output_file_made = False
contaminat_output_file = sys.argv[4]
contaminat_output_seqs = []
contaminat_output_file_made = False

start_all = clock()
with open(BLAT_tab_file, "r") as tabfile:
    contaminated_seqs = []
    query_seq = ""
    for line in tabfile:
        if len(line) < 2:
            continue
        line_parts = line.split("\t")
        if query_seq == line_parts[0]:
            continue
        else:
            query_seq = line_parts[0]
            contaminated_seqs.append(query_seq)
            
    for seq in input_seqs:
        if seq not in contaminated_seqs:
            output_seqs.append(input_seqs[seq])
            if len(output_seqs) > 100000:
                if not output_file_made:
                    output_file_made = True
                    with open(output_file, "w") as outfile:
                        SeqIO.write(output_seqs, outfile, "fastq")
                else:
                    with open(output_file, "a") as outfile:
                        SeqIO.write(output_seqs, outfile, "fastq")
                output_seqs = []
        else:
            contaminat_output_seqs.append(input_seqs[seq])
            if len(contaminat_output_seqs) > 100000:
                if not contaminat_output_file_made:
                    contaminat_output_file_made = True
                    with open(contaminat_output_file, "w") as outfile:
                        SeqIO.write(contaminat_output_seqs, outfile, "fastq")
                else:
                    with open(contaminat_output_file, "a") as outfile:
                        SeqIO.write(contaminat_output_seqs, outfile, "fastq")
                contaminat_output_seqs = []
                
                
    if len(output_seqs) > 0:
        if not output_file_made:
            output_file_made = True
            with open(output_file, "w") as outfile:
                SeqIO.write(output_seqs, outfile, "fastq")
        else:
            with open(output_file, "a") as outfile:
                SeqIO.write(output_seqs, outfile, "fastq")
    if len(contaminat_output_seqs) > 0:
        if not contaminat_output_file_made:
            contaminat_output_file_made = True
            with open(contaminat_output_file, "w") as outfile:
                SeqIO.write(contaminat_output_seqs, outfile, "fastq")
        else:
            with open(contaminat_output_file, "a") as outfile:
                SeqIO.write(contaminat_output_seqs, outfile, "fastq")
                
end_all = clock()
#print "--------------------------------------"
#print "BLAT containment filter"
#print "==============================="
#print "total run:", end_all - start_all, "s"          