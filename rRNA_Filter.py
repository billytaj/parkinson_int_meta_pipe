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

Barrnap = "/home/j/jparkins/mobolaji/Tools/Barrnap/bin/barrnap"
Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch"
Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"

sequence_file = sys.argv[1]
sequence_file_fasta = os.path.splitext(sequence_file)[0] + ".fasta"
mRNA_file = sys.argv[2]
rRNA_file = sys.argv[3]

sequences = list(SeqIO.parse(sequence_file, "fastq"))
with open(sequence_file_fasta, "w") as out:
    SeqIO.write(sequences, out, "fasta")

Barrnap_rRNA_IDs = set()
Barrnap_out = os.path.splitext(sequence_file)[0] + "_n_rRNA_barrnap.fasta"
Infernal_rRNA_IDs = set()
Infernal_out = os.path.splitext(sequence_file)[0] + "_infernal.tab"

mRNA_seqs = set()
rRNA_seqs = set()

arc = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "arc", "--reject", "0.01", sequence_file_fasta])
bac = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "bac", "--reject", "0.01", sequence_file_fasta])
euk = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "euk", "--reject", "0.01", sequence_file_fasta])
mito = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "mito", "--reject", "0.01", sequence_file_fasta])

for line in arc.splitlines():
    if not line.startswith("#"):
        Barrnap_rRNA_IDs.add(line.split("\t")[0])

for line in bac.splitlines():
    if not line.startswith("#"):
        Barrnap_rRNA_IDs.add(line.split("\t")[0])

for line in euk.splitlines():
    if not line.startswith("#"):
        Barrnap_rRNA_IDs.add(line.split("\t")[0])

for line in mito.splitlines():
    if not line.startswith("#"):
        Barrnap_rRNA_IDs.add(line.split("\t")[0])

for sequence in sequences:
    if sequence.id in Barrnap_rRNA_IDs:
        rRNA_seqs.add(sequence)
    else:
        mRNA_seqs.add(sequence)

with open(Barrnap_out, "w") as out:
    SeqIO.write(list(mRNA_seqs), out, "fasta")

subprocess.call([Infernal, "-o", "/dev/null", "--tblout", Infernal_out, "--anytrunc", "--rfam", "-E", "0.001", Rfam, Barrnap_out])

with open(Infernal_out, "r") as infile_read:
    for line in infile_read:
        if not line.startswith("#") and len(line) > 10:
            Infernal_rRNA_IDs.add(line[:line.find(" ")].strip())

mRNA_seqs = set()

for sequence in sequences:
    if sequence.id in Infernal_rRNA_IDs or sequence.id in Barrnap_rRNA_IDs:
        rRNA_seqs.add(sequence)
    else:
        mRNA_seqs.add(sequence)

with open(mRNA_file, "w") as out:
    SeqIO.write(list(mRNA_seqs), out, "fastq")

with open(rRNA_file, "w") as out:
    SeqIO.write(list(rRNA_seqs), out, "fastq")

if len(sys.argv) > 4:
    sequence_file_pair = sys.argv[4]
    sequences_pair = list(SeqIO.parse(sequence_file_pair, "fastq"))
    mRNA_file_pair = sys.argv[5]
    rRNA_file_pair = sys.argv[6]
    mRNA_seqs = set()
    rRNA_seqs = set()
    for sequence in sequences_pair:
        if sequence.id in Infernal_rRNA_IDs or sequence.id in Barrnap_rRNA_IDs:
            rRNA_seqs.add(sequence)
        else:
            mRNA_seqs.add(sequence)

    with open(mRNA_file_pair, "w") as out:
        SeqIO.write(list(mRNA_seqs), out, "fastq")

    with open(rRNA_file_pair, "w") as out:
        SeqIO.write(list(rRNA_seqs), out, "fastq")

### ADD PERCENT ID FILTER 90%

end_all = clock()
print("rRNA filter")
print("============================================")
print("total runtime", end_all - start_all, "s")