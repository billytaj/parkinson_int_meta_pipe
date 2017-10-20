#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re



def sortbyscore(line):
    return line[11]
	
	
def gene_map(tsv, unmapped, read_seqs, gene2read_map, contig2read_map, mapped_reads):
        with open(tsv, "r") as tabfile:
            query = ""
            identity_cutoff = 85
            length_cutoff = 0.65
            score_cutoff = 60
            Hits = []
			b_contig_hit = False
            for line in tabfile:
                if len(line) < 2:
                    continue
                else:
                    Hits.append(line.split("\t"))
            Sorted_Hits = sorted(Hits, key = sortbyscore)
            for line in Sorted_Hits:
                if query in contig2read_map:
                    b_contig_hit = True
                else:
                    b_contig_hit = False
                if query == line[0]:
                    continue
                else:
                    query = line[0]
                    db_match = line[1]
                    seq_identity = line[2]
                    align_len = line[3]
                    score = line[11]
                if float(seq_identity) > float(identity_cutoff):
                    if align_len > len(read_seqs[query].seq) * length_cutoff:
                        if float(score) > float(score_cutoff):
                            if db_match in gene2read_map:
                                if b_contig_hit:
                                    for read in contig2read_map[query]:
                                        if read not in gene2read_map[db_match]:
                                            if read not in mapped_reads:
                                                gene2read_map[db_match].append(read)
                                                mapped_reads.add(read)
                                else:
                                    if query not in gene2read_map[db_match]:
                                        gene2read_map[db_match].append(query)
                                        mapped_reads.add(query)
                            else:
                                if b_contig_hit:
                                    read_count = 0
                                    for read in contig2read_map[query]:
                                        if read not in mapped_reads:
                                            mapped_reads.add(read)
                                            read_count += 1
                                            if read_count == 1:
                                                gene2read_map[db_match] = [read]
                                            elif read_count > 1:
                                                gene2read_map[db_match].append(read)
                                else:
                                    gene2read_map[db_match] = [query]
                                    mapped_reads.add(query)
                            continue
                unmapped.add(query)	

				
def main():		
	
	if(len(sys.argv) == 2):
		if(sys.argv[1] == "arg_list"):
			arg_list  = []
			gene2read_map = {}
			mapped_reads = set()
			unmapped_reads = set()
			unmapped_seqs = []
			reads_count = 0
			genes = []
			
			with open (sys.argv[1], 'r') as arg_file:
				arg_list = arg_file.readlines()
			arg_file.closed
			
			dna_db_path = arg_list[0]
			contig_to_read_path = arg_list[1]
			gene_file_to_read = arg_list[2]			
			gene_map = arg_list[3]
			contigs_not_BWA_path = arg_list[4]
			output_path = arg_list[5]
			contig_not_aligned = arg_list[6]
			
			read_seqs = SeqIO.index(contigs_not_BWA_path, os.path.splitext(contigs_not_BWA_path)[1][1:])
			
			with open(contig_to_read_path, "r") as mapping:
			for line in mapping:
				if len(line) > 5:
					entry = line.strip("\n").split("\t")
					contig2read_map[entry[0]] = entry[2:]

			gene2read_map = {}
			mapped_reads = set()

			with open(gene2read_file, "r") as mapping:
				for line in mapping:
					if len(line) > 5:
						entry = line.split("\t")
						gene2read_map[entry[0]] = entry[3:]
			
			
			gene_map(BLAT_tab_file, unmapped_reads)
			
			
			# for each seq in the not_BWA section...
			for read in read_seqs:
				if read not in unmapped_reads:
					for gene in gene2read_map:
						if read in gene2read_map[gene]:
							break
					else:
						unmapped_reads.add(read)

			for read in unmapped_reads:
				unmapped_seqs.append(read_seqs[read])
			with open(output_file, "w") as outfile:
				SeqIO.write(unmapped_seqs, outfile, "fasta")
				
			with open(gene_file_to_read, "w") as out_map:
				for record in SeqIO.parse(DNA_DB, "fasta"):
					if record.id in gene2read_map:
						genes.append(record)
						out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
						for read in gene2read_map[record.id]:
							out_map.write("\t" + read.strip("\n"))
							reads_count += 1
						else:
							out_map.write("\n")
			with open(gene_file, "w") as outfile:
				SeqIO.write(genes, outfile, "fasta")	
			
		else:
			print("not accepted")
	
	elif(len(sys.argv) > 2):
		#old way of doing stuff
		DNA_DB = sys.argv[1]
		contig2read_file = sys.argv[2]
		gene2read_file = sys.argv[3]
		gene_file = sys.argv[4]

		contig2read_map = {}

		with open(contig2read_file, "r") as mapping:
			for line in mapping:
				if len(line) > 5:
					entry = line.strip("\n").split("\t")
					contig2read_map[entry[0]] = entry[2:]

		gene2read_map = {}
		mapped_reads = set()

		with open(gene2read_file, "r") as mapping:
			for line in mapping:
				if len(line) > 5:
					entry = line.split("\t")
					gene2read_map[entry[0]] = entry[3:]
		
		# for each grouping
		for x in range((len(sys.argv) - 5) / 3):
			read_file = sys.argv[3 * x + 5]
			read_seqs = SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
			BLAT_tab_file = sys.argv[3 * x + 6]
			output_file = sys.argv[3 * x + 7]

			unmapped_reads = set()
			unmapped_seqs = []

			

			gene_map(BLAT_tab_file, unmapped_reads)
			# for each seq in the not_BWA section...
			for read in read_seqs:
				if read not in unmapped_reads:
					for gene in gene2read_map:
						if read in gene2read_map[gene]:
							break
					else:
						unmapped_reads.add(read)

			for read in unmapped_reads:
				unmapped_seqs.append(read_seqs[read])
			with open(output_file, "w") as outfile:
				SeqIO.write(unmapped_seqs, outfile, "fasta")

		reads_count = 0
		genes = []
		with open(gene2read_file, "w") as out_map:
			for record in SeqIO.parse(DNA_DB, "fasta"):
				if record.id in gene2read_map:
					genes.append(record)
					out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
					for read in gene2read_map[record.id]:
						out_map.write("\t" + read.strip("\n"))
						reads_count += 1
					else:
						out_map.write("\n")
		with open(gene_file, "w") as outfile:
			SeqIO.write(genes, outfile, "fasta")

		print (str(reads_count) + " reads were mapped with BWA and BLAT")
		print ("Reads mapped to " + str(len(genes)) + " genes.")
		
	else:
		print("not enough args")


if __name__ == "__main__":
	main()
	