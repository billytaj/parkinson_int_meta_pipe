#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock

start_all = clock()

Read1_file = sys.argv[1]
Read1_seqs = SeqIO.index(Read1_file, "fastq")
Read1_out = sys.argv[2]
Read2_file = sys.argv[3]
Read2_seqs = SeqIO.index(Read2_file, "fastq")
Read2_out = sys.argv[4]
Unpaired_file = sys.argv[5]

Paired1_seqs = []
Paired2_seqs = []
Unpaired_seqs = []

file_made_paired = False

for seq in Read1_seqs:
    if seq in Read2_seqs:
        Paired1_seqs.append(Read1_seqs[seq])
        Paired2_seqs.append(Read2_seqs[seq])
        if len(Paired1_seqs) > 100000:
            if not file_made_paired:
                file_made_paired = True
                with open(Read1_out, "w") as out:
                    SeqIO.write(Paired1_seqs, out, "fastq")
                with open(Read2_out, "w") as out:
                    SeqIO.write(Paired2_seqs, out, "fastq")
                Paired1_seqs = []
                Paired2_seqs = []
            else:
                with open(Read1_out, "a") as out:
                    SeqIO.write(Paired1_seqs, out, "fastq")
                with open(Read2_out, "a") as out:
                    SeqIO.write(Paired2_seqs, out, "fastq")
                Paired1_seqs = []
                Paired2_seqs = []
    else:
        Unpaired_seqs.append(Read1_seqs[seq])
        if len(Unpaired_seqs) > 100000:
            with open(Unpaired_file, "a") as out:
                SeqIO.write(Unpaired_seqs, out, "fastq")
            Unpaired_seqs = []
if len(Paired1_seqs) > 0:
    if not file_made_paired:
        file_made_paired = True
        with open(Read1_out, "w") as out:
            SeqIO.write(Paired1_seqs, out, "fastq")
        with open(Read2_out, "w") as out:
            SeqIO.write(Paired2_seqs, out, "fastq")
    else:
        with open(Read1_out, "a") as out:
            SeqIO.write(Paired1_seqs, out, "fastq")
        with open(Read2_out, "a") as out:
            SeqIO.write(Paired2_seqs, out, "fastq")

for seq in Read2_seqs:
    if seq not in Read1_seqs:
        Unpaired_seqs.append(Read2_seqs[seq])
        if len(Unpaired_seqs) > 100000:
            with open(Unpaired_file, "a") as out:
                SeqIO.write(Unpaired_seqs, out, "fastq")
            Unpaired_seqs = []
if len(Unpaired_seqs) > 0:
    with open(Unpaired_file, "a") as out:
        SeqIO.write(Unpaired_seqs, out, "fastq")

end_all = clock()
print("Paired reads filter")
print("===========================================")
print("total run time:", end_all - start_all, "s")

