#@@ -1,258 +0,0 @@
#!/usr/bin/env python 
 
import sys 
import os 
import os.path 
import shutil 
import subprocess 
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord 
import re 
import time 
from datetime import datetime as dt 
 
def sortbyscore(line): 
    return line[11] 
     
     
def gene_map(blat_results, unmapped, read_seqs, gene2read_map, contig2read_map, mapped_reads): 
    start_internal_gene_map_time = time.clock() 
     
    align_len = 0 
    seq_identity = 0 
    score = 0 
    db_match = "" 
 
    with open(blat_results, "r") as tabfile: 
        query = "" 
        identity_cutoff = 85.0 
        length_cutoff = 0.65 
        score_cutoff = 60.0 
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
                seq_identity = float(line[2]) 
                align_len = float(line[3]) 
                score = float(line[11]) 
                #print("sorted hit out:", line) 
                #sys.exit() 
            if ((seq_identity > identity_cutoff)  
            and (align_len > (len(read_seqs[query].seq) * length_cutoff))  
            and (score > score_cutoff)): 
                if db_match in gene2read_map: 
                    if b_contig_hit: 
                        for read in contig2read_map[query]: 
                            mapped_reads.add(read) 
                            if (read not in gene2read_map[db_match]): 
                                gene2read_map[db_match].append(read) 
                                 
                    else: 
                        if query not in gene2read_map[db_match]: 
                            gene2read_map[db_match].append(query) 
                            mapped_reads.add(query) 
                else: 
                    if b_contig_hit: 
                        read_count = 0 
                        for read in contig2read_map[query]: 
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
    end_internal_gene_map_time = time.clock() 
    print("internal gene map time:", end_internal_gene_map_time - start_internal_gene_map_time, "s") 
                 
def main():      
     
    if(len(sys.argv) == 3): 
        if(sys.argv[1] == "arg_list"): 
            arg_list  = [] 
            gene2read_map = {} 
            mapped_reads = set() 
            unmapped_contig_reads = set() 
            unmapped_contig_seqs = [] 
             
            unmapped_unpaired_reads = set() 
            unmapped_unpaired_seqs = [] 
             
            unmapped_pair_1_reads = set() 
            unmapped_pair_1_seqs = [] 
             
            unmapped_pair_2_reads = set() 
            unmapped_pair_2_seqs = [] 
            reads_count = 0 
            genes = [] 
            contig2read_map = {} 
 
            with open (sys.argv[2], 'r') as arg_file: 
                arg_list = arg_file.readlines() 
            arg_file.closed 
             
            dna_db_path = arg_list[0].split('\n')[0] 
            contig_to_read_path = arg_list[1].split('\n')[0] 
            gene_file_to_read = arg_list[2]         .split('\n')[0] 
            final_gene_map = arg_list[3].split('\n')[0] # just gets dumped-to.  doesn't matter what's in it, really 
            #---------------------------- 
            #repeats here in batches of 2 
            contigs_not_BWA_path = arg_list[4].split('\n')[0] 
            contig_blat_results = arg_list[5].split('\n')[0] 
             
            unpaired_unmapped_not_BWA_path = arg_list[6].split('\n')[0] 
            unpaired_unmapped_not_BWA_blat_results = arg_list[7].split('\n')[0] 
             
            pair_1_not_BWA_path = arg_list[8].split('\n')[0] 
            pair_1_not_BWA_blat_results = arg_list[9].split('\n')[0] 
             
            pair_2_not_BWA_path = arg_list[10].split('\n')[0] 
            pair_2_not_BWA_blat_results = arg_list[11].split('\n')[0] 
             
             
            #------------------------------------------------ 
            # open the files and put it through the gene map 
             
            print("contigs not BWA loc:", contigs_not_BWA_path) 
            print("format:", os.path.splitext(contigs_not_BWA_path)[1][1:]) 
            #read_seqs = SeqIO.index(contigs_not_BWA_path, os.path.splitext(contigs_not_BWA_path)[1][1:]) 
            contig_read_seqs = SeqIO.index(contigs_not_BWA_path, 'fasta') 
             
            unpaired_read_seqs = SeqIO.index(unpaired_unmapped_not_BWA_path, 'fasta') 
            pair_1_read_seqs = SeqIO.index(pair_1_not_BWA_path, 'fasta') 
            pair_2_read_seqs = SeqIO.index(pair_2_not_BWA_path, 'fasta') 
             
            start_read_contig_time = time.clock() 
            with open(contig_to_read_path, "r") as contig_mapping: 
                for line in contig_mapping: 
                    if len(line) > 5: 
                        entry = line.strip("\n").split("\t") 
                        contig2read_map[entry[0]] = entry[2:] 
 
            end_read_contig_time = time.clock() 
            start_read_gene_map = time.clock() 
            with open(gene_file_to_read, "r") as gene_mapping: 
                for line in gene_mapping: 
                    if len(line) > 5: 
                        entry = line.split("\t") 
                        gene2read_map[entry[0]] = entry[3:] 
            end_read_gene_map = time.clock() 
           
            start_gene_map_time = time.clock() 
            gene_map(contig_blat_results, unmapped_contig_reads, contig_read_seqs, gene2read_map, contig2read_map, mapped_reads) 
            gene_map(unpaired_unmapped_not_BWA_blat_results, unmapped_unpaired_reads, unpaired_read_seqs, gene2read_map, contig2read_map, mapped_reads) 
            gene_map(pair_1_not_BWA_blat_results, unmapped_pair_1_reads, pair_1_read_seqs, gene2read_map, contig2read_map, mapped_reads) 
            gene_map(pair_2_not_BWA_blat_results, unmapped_pair_2_reads, pair_2_read_seqs, gene2read_map, contig2read_map, mapped_reads) 
             
            end_gene_map_time = time.clock() 
             
            # for each seq in the not_BWA section... 
            start_first_write_time = time.clock() 
            for read in contig_read_seqs: 
                if(read not in unmapped_contig_reads): 
                    #unmapped_reads is empty at start.  Only appends on the first hit of a new value 
                    unmapped_contig_seqs.append(contig_read_seqs[read]) 
                unmapped_contig_reads.add(read) 
             
            for read in unpaired_read_seqs: 
                if(read not in unmapped_unpaired_reads): 
                    #unmapped_reads is empty at start.  Only appends on the first hit of a new value 
                    unmapped_unpaired_seqs.append(unpaired_read_seqs[read]) 
                unmapped_unpaired_reads.add(read) 
             
            for read in pair_1_read_seqs: 
                if(read not in unmapped_pair_1_reads): 
                    #unmapped_reads is empty at start.  Only appends on the first hit of a new value 
                    unmapped_pair_1_seqs.append(pair_1_read_seqs[read]) 
                unmapped_pair_1_reads.add(read) 
             
            for read in pair_2_read_seqs: 
                if(read not in unmapped_pair_2_reads): 
                    #unmapped_reads is empty at start.  Only appends on the first hit of a new value 
                    unmapped_pair_2_seqs.append(pair_2_read_seqs[read]) 
                unmapped_pair_2_reads.add(read) 
                 
             
            with open("contigs_not_aligned.fasta", "w+") as outfile: 
                SeqIO.write(unmapped_contig_seqs, outfile, "fasta") 
             
            with open("unpaired_not_aligned.fasta", "w+") as outfile: 
                SeqIO.write(unmapped_unpaired_seqs, outfile, "fasta") 
                 
            with open("pair_1_not_aligned.fasta", "w+") as outfile: 
                SeqIO.write(unmapped_pair_1_seqs, outfile, "fasta") 
 
            with open("pair_2_not_aligned.fasta", "w+") as outfile: 
                SeqIO.write(unmapped_pair_2_seqs, outfile, "fasta") 
             
            end_first_write_time = time.clock()     
                 
            #-------------------------------------------- 
             
 
 
            #----------------------------------------------------------------- 
            #this stuff doesn't get repeated 
            start_write_time = time.clock()   
            with open(gene_file_to_read, "w") as out_map: 
                diagnostic_len_flag = True 
                dna_file = SeqIO.parse(dna_db_path, "fasta") 
                dna_record_count = 0 
                gene_map_record_count = 0 
                for record in dna_file: 
                    dna_record_count += 1 
                    #print("record id:", record.id) 
                    if record.id in gene2read_map: 
                        genes.append(record) 
                        out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id]))) 
                        for read in gene2read_map[record.id]: 
                            out_map.write("\t" + read.strip("\n")) 
                            gene_map_record_count += 1 
                    else: 
                        out_map.write("\n") 
            end_write_time = time.clock()  
            final_name = final_gene_map.split(".")[0] + "_" + dna_db_path.split(".")[0] + ".fna" 
            with open(final_name, "w+") as outfile: 
                SeqIO.write(genes, outfile, "fasta")     
            print("total number of dna in dna_db:", dna_record_count) 
            print("total number of gene map records to scan through:", gene_map_record_count) 
            print("PROGRAM DONE") 
            print("output file:", final_name) 
            print("read contig time:", end_read_contig_time - start_read_contig_time, "s") 
            print("final gene map time", end_gene_map_time - start_gene_map_time, "seconds") 
            print("contig not aligned write time:", end_first_write_time - start_first_write_time, "s") 
            print("gene map write time:", end_write_time - start_write_time, "s") 
             
        else: 
            print("not accepted") 
    else: 
        print("not enough args") 
 
 
if __name__ == "__main__": 
    print("job started at:", dt.today())  
    start_main_time = time.clock() 
    main() 
    end_main_time = time.clock() 
    print("total program time:", end_main_time - start_main_time, "seconds")