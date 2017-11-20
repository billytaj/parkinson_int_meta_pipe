#!/usr/bin/env python

import os
import os.path
import sys
import time
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Input_Dir = sys.argv[1]
Split_Dir = sys.argv[2]

try:
    os.mkdir(Split_Dir)
except:
    shutil.rmtree(Split_Dir)
    os.mkdir(Split_Dir)

split_len = 1000

Nucleotides = set(["A", "T", "G", "C", "-", "N"])

def check_format(seq_obj):
	seq_format = ""
	for sequence in seq_obj:
		if seq_format == "RNA" or seq_format == "PROT":
			break
		for character in sequence.seq:
			if character.upper() in Nucleotides:
				continue
			elif character.upper() == "U":
				seq_format = "RNA"
				break
			else:
				seq_format = "PROT"
				break
		else:
			seq_format = "DNA"
			break
	return seq_format

for file in os.listdir(Input_Dir):
	if os.path.isdir(os.path.join(Input_Dir, file)):
		continue
	sequences = list(SeqIO.parse(os.path.join(Input_Dir, file), "fasta"))
	file_lenght = len(sequences)
	seq_format = check_format(sequences)
	if seq_format == "DNA":
		for index, sequence in enumerate(sequences):
			try:
				sequences[index] = SeqRecord(seq = sequence.seq.translate(stop_symbol=""), id=sequence.id, description=sequence.description)
			except:
				pass
		else:
			seq_format = "PROT"
	elif seq_format == "RNA":
		print "This is an mRNA sequence!!!"
	if file_lenght > split_len:
		for i in range(file_lenght/split_len):
			file_count = i
			start = file_count*split_len
			file_count += 1
			stop = file_count*split_len
			temp_sequences = []
			for index, sequence in enumerate(sequences):
				if index >= start and index < stop:
					temp_sequences.append(sequence)
				elif index == stop:
					break
			with open(os.path.join(Split_Dir, os.path.splitext(file)[0] + "_split_" + str(file_count) + ".faa"), "w") as out_faa:
				SeqIO.write(temp_sequences, out_faa, "fasta")
	else:
		with open(os.path.join(Split_Dir, os.path.splitext(file)[0] + "_split_0.faa"), "w") as out_faa:
			SeqIO.write(sequences, out_faa, "fasta")