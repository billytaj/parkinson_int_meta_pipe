#!/usr/bin/env python
# this script calls infernal.
# this script gets called by a higher power
# Most of the work here seems to just want to make Barrnap play nice.
# But after further discussion, it's just a lightweight, speedy mini-filter to take the load off of infernal (which was supposed to be slow)
import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock
import mt_pipe_paths as mpp

start_all = clock()



sequence_file = sys.argv[1]#in
sequence_file_fasta = os.path.splitext(sequence_file)[0] + ".fasta"
mRNA_file = sys.argv[2] #out
rRNA_file = sys.argv[3] #out

start_write_sequence_fasta = clock()
sequences = list(SeqIO.parse(sequence_file, "fastq"))
with open(sequence_file_fasta, "w") as out:
    SeqIO.write(sequences, out, "fasta")
end_write_sequence_fasta = clock()
    
Barrnap_rRNA_IDs = set()
Barrnap_out = os.path.splitext(sequence_file)[0] + "_n_rRNA_barrnap.fasta"
Infernal_rRNA_IDs = set()
Infernal_out = os.path.splitext(sequence_file)[0] + "_infernal.tab"

mRNA_seqs = set()
rRNA_seqs = set()
start_barrnap_arc = clock()
arc = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "arc", "--reject", "0.01", sequence_file_fasta])
end_barrnap_arc = clock()
start_barrnap_bac = clock()
bac = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "bac", "--reject", "0.01", sequence_file_fasta])
end_barrnap_bac = clock()
start_barrnap_euk = clock()
euk = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "euk", "--reject", "0.01", sequence_file_fasta])
end_barrnap_euk = clock()
start_barrnap_mito = clock()
mito = subprocess.check_output([Barrnap, "--quiet", "--kingdom", "mito", "--reject", "0.01", sequence_file_fasta])
end_barrnap_mito = clock()

start_barrnap_text_manip = clock()
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
end_barrnap_text_manip = clock()

start_write_mrna_fasta = clock()        
with open(Barrnap_out, "w") as out:
    SeqIO.write(list(mRNA_seqs), out, "fasta")
end_write_mrna_fasta = clock()
    
start_infernal = clock()
subprocess.call([Infernal, "-o", "/dev/null", "--tblout", Infernal_out, "--anytrunc", "--rfam", "-E", "0.001", Rfam, Barrnap_out])
end_infernal = clock()

start_infernal_result_parse = clock()
with open(Infernal_out, "r") as infile_read:
    for line in infile_read:
        if not line.startswith("#") and len(line) > 10:
            Infernal_rRNA_IDs.add(line[:line.find(" ")].strip())
end_infernal_result_parse = clock()
            
mRNA_seqs = set()

for sequence in sequences:
    if sequence.id in Infernal_rRNA_IDs or sequence.id in Barrnap_rRNA_IDs:
        rRNA_seqs.add(sequence)
    else:
        mRNA_seqs.add(sequence)

start_write_mrna_fastq = clock()        
with open(mRNA_file, "w") as out:
    SeqIO.write(list(mRNA_seqs), out, "fastq")
end_write_mrna_fastq = clock()
start_write_rrna_fastq = clock()
with open(rRNA_file, "w") as out:
    SeqIO.write(list(rRNA_seqs), out, "fastq")
end_write_rrna_fastq = clock()
    
start_alternate_write_fastq = clock()
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
end_alternate_write_fastq = clock()
### ADD PERCENT ID FILTER 90%

end_all = clock()

with open(os.path.splitext(sequence_file)[0] + "_rna_filter_profile.txt", 'w+') as profile:
    profile.write("rRNA filter\n")
    profile.write("============================================\n")
    profile.write("total runtime: " +  str(end_all - start_all) + "s\n")
    profile.write("sequence file write time: " + str(end_write_sequence_fasta - start_write_sequence_fasta) + "s\n")
    profile.write("barrnap arc time: " + str(end_barrnap_arc - start_barrnap_arc) + "s\n")
    profile.write("barrnap bac time: " + str(end_barrnap_bac - start_barrnap_bac) + "s\n")
    profile.write("barrnap euk time: " + str(end_barrnap_euk - start_barrnap_euk) + "s\n")
    profile.write("barrnap mito time: " + str(end_barrnap_mito - start_barrnap_mito) + "s\n")
    profile.write("barrnap text manipulation: " + str(end_barrnap_text_manip - start_barrnap_text_manip) + "s\n")
    profile.write("writing mrna FASTA: " + str(end_write_mrna_fasta - start_write_mrna_fasta) + "s\n")
    profile.write("INFERNAL TIME: " + str(end_infernal - start_infernal) + "s\n")
    profile.write("parsing the INFERNAL output: " + str(end_infernal_result_parse - start_infernal_result_parse) + "s\n")
    profile.write("writing mrna FASTQ: " + str(end_write_mrna_fastq - start_write_mrna_fastq) + "s\n")
    profile.write("writing rrna FASTQ: " + str(end_write_rrna_fastq - start_write_rrna_fastq) + "s\n")
    profile.write("writing alternate FASTQs: " + str(end_alternate_write_fastq - start_alternate_write_fastq) + "s\n")
    profile.close()


