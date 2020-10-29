#!/usr/bin/env python
#This version of the BLAT PP code is only concerned with 1 set of data.
#it will take in a BLAT output, and seek to form the reads into a map of genes to constituent reads
#it will not import the old gene map.  Merging will be done at the end of the GA phase.  
#
#oct 22, 2020:
#New in v2: top-hit! 

import os
import os.path
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from datetime import datetime as dt
import multiprocessing as mp
from shutil import copyfile
import time
from math import ceil
def import_contig_map(contig2read_file):
    # make dict of contigID<->readsID(s):
    contig2read_map= {}
    with open(contig2read_file,"r") as mapping:
        for line in mapping:
            if len(line)>5:                             # line starts with 'NODE_'
                entry= line.strip("\n").split("\t")     # break tab-separated into list
                contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
    return contig2read_map

def extract_true_read_length(read):
    real_length = (len(read) - 10)/2
    if(real_length % 1 != 0):
        sys.exit("reads aren't of even length")
    return real_length

def import_gene_map(gene2read_file):
    # make dict of BWA-aligned geneID<->readID(s):
    #BWAreads = []                                        # DEBUG
    gene2read_map = defaultdict(list)                    # Dict of BWA&BLAT-aligned geneID<->readID(s)
    mapped_reads = set()
    with open(gene2read_file,"r") as mapping:           #  initialized w BWA-alignments.
        for line in mapping:
            if len(line)>5:                             # line at least 5 characeters?
                entry = line.strip("\n").split("\t")
                gene2read_map[entry[0]] = entry[3:]      # key=geneID, value=list of readID(s)
                                                        # (using [:] syntax ensures a list, even if one ele)
                mapped_reads.update(entry[3:])
                #BWAreads.extend(entry[3:])              # DEBUG
    return gene2read_map, mapped_reads#, BWAreads
    

#####################################
# FUNCTIONS:

def check_file_safety(file_name):
    if(os.path.exists(file_name)):
        if(os.path.getsize(file_name) == 0):
            #copyfile(gene2read_file, new_gene2read_file)
            print(file_name, "is unsafe.  skipping")
            return False
            #sys.exit(DMD_tab_file_1 + " -> DMD tab file 1 is empty.  aborting")
        else:
            print(file_name, "exists and is safe")
            return True
    else:
        #sys.exit(DMD_tab_file_1 + " -> DMD tab file 1 is missing.  aborting")
        print(file_name, "is missing. skipping")
        return False


# read .blatout file and acquire read length:
def get_blat_details(blat_in, reads_in):                          # List of lists [blatout info, read length]=
                                                        #  read_aligned(.blatout file, dict contig/readID<->SeqRecord)
    
    seqrec = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:])
    
    #for item in seqrec:
    #    print(item, seqrec[item])
        
        
        
    # get info from .blatout file:
    print ('Reading ' + str(os.path.basename(blat_in)) + '.')
    with open(blat_in,"r") as tabfile:
        hits = []                                        # List of lists containing .blatout fields.
        for line in tabfile:                            # In the .blatout file:
            if len(line)>=2:                            # If length of line >= 2,
                info_list= line.strip("\n").split("\t") #  make a list from .blatout tab-delimited fields,
                query= info_list[0]                     #  get the queryID= contig/readID,
                info_list.append(len(seqrec[query].seq))#  add query length (int) to end of list,
                
                #info_list.append(len(seqrec[query]))#  add query length (int) to end of list,
                hits.append(info_list)                  #  and append the list to the hits list.

    # return info:
    print(dt.today(), "hits in blatout:", len(hits))
    
    #send back the hits based on e-value
    hits = sorted(hits, key=sortbyscore, reverse = True)
    
    return hits, seqrec




# sort by score:
# sort by the bitscore. bitscore > e_Value.  bitscores don't change on database sizes
def sortbyscore(line):
    return float(line[11])

# add BLAT-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:

def make_gene_map(hits, contig2read_map):                         # fail-mapped contig/readIDs=
    # sort alignment list by high score:
    gene2read_map = dict()
    mapped = set()
    unmapped = set()                         # Set of unmapped contig/readIDs.
    query_details_dict = dict()

    # BLAT threshold:
    identity_cutoff= 85
    #identity_cutoff= 80
    length_cutoff= 0.65
    true_length_cutoff = 0.85
    score_cutoff= 60

    repeat_reads = 0
    disagreements = 0
    one_sided_alignments = 0
    rejected = 0
    accepted = 0
    total_reads = 0
    rejected_score = 0
    rejected_align = 0
    rejected_id = 0
    # loop through sorted BLAT-aligned reads/contigs:
    for line in hits:
        total_reads += 1
        inner_details_dict = dict()
        # "ON-THE-FLY" filter:
        
        # extract & store data:
        query = line[0]                      # queryID= contig/readID
        db_match = line[1]                   # geneID
        seq_identity = float(line[2])        # sequence identity
        align_len = int(line[3])             # alignment length
        score = float(line[11])              # score
        seq_len = int(line[12])              # query sequence length
        bitscore = float(line[11])
        
        true_seq_len = ceil((seq_len - 10)/2)
        if(true_seq_len % 1 != 0):  
            print(dt.today(), "reads not of equal length", query, seq_len, true_seq_len)
            sys.exit("death")
        
        # does query alignment meet threshold?:
        # (identity, length, score):
        if seq_identity<=identity_cutoff  or score<=score_cutoff or align_len<=true_seq_len*true_length_cutoff:
            #unmapped.add(query)             # if threshold is failed, add to unmapped set and
            if(seq_identity <= identity_cutoff):
                rejected_id += 1
            elif(align_len<=true_seq_len*true_length_cutoff):
                rejected_align += 1
                if(query not in query_details_dict):
                    print(dt.today(), query, "rejected align:", align_len/true_seq_len)
                time.sleep(0.05)
            elif(score <= score_cutoff):
                rejected_score += 1
                

            inner_details_dict = dict()
            inner_details_dict["match"] = False
            query_details_dict[query] = inner_details_dict
            rejected += 1
            

        else:
            
            # is query a contig?:
            if query in contig2read_map:        # If query is found in contig list,
                contig= True                    #  then mark as contig,
            else:                               #  if not,
                contig= False                   #  mark as not a contig.
                
                
            # is this alignment the highest-score match for that query?
            if query in query_details_dict: #query2gene_map:         # If alignment previously found for this contig/read,
                repeat_reads += 1
                inner_details_dict = query_details_dict[query]
                if(inner_details_dict["match"]):
                    old_bitscore = inner_details_dict["bitscore"]
                    disagreements += 1
                    if(old_bitscore < bitscore):
                        print("old:", old_bitscore, "new:", bitscore)
                        #update, but this shouldn't happen
                        inner_details_dict["gene"] = db_match
                        inner_details_dict["bitscore"] = bitscore
                        inner_details_dict["is_contig"] = contig
                else:
                    rejected -= 1
                    accepted += 1
                    one_sided_alignments += 1
                    inner_details_dict["match"] = True
                    inner_details_dict["is_contig"] = contig
                    inner_details_dict["gene"] = db_match
                    inner_details_dict["bitscore"] = bitscore
                
                                            
                              # skip to the next query.
            else:
                #new hit
                accepted += 1
                inner_details_dict = dict()
                inner_details_dict["match"] = True
                inner_details_dict["is_contig"] = contig
                inner_details_dict["gene"] = db_match
                inner_details_dict["bitscore"] = bitscore
                query_details_dict[query] = inner_details_dict
            #queryIScontig[query] = contig        # Store contig (T/F) info.
            #query2gene_map[query].add(db_match) # Collect all aligned genes for contig/read;
                                                #  query2gene_map[query] is a set, although there should be
                                                #  no more than one gene alignement getting through filter.
            
    # EXPAND HERE to deal with multiple high-score genes for a read.
    print(dt.today(), "repeat reads:", repeat_reads)
    print(dt.today(), "disagreements:", disagreements)
    print(dt.today(), "one-sided alignments:", one_sided_alignments)
    print(dt.today(), "accepted:", accepted, "/", total_reads)
    print(dt.today(), "rejected:", rejected, "/", total_reads)
    print(dt.today(), "rejected by seq id:", rejected_id, "/", total_reads)
    print(dt.today(), "rejected by align:", rejected_align, "/", total_reads)
    print(dt.today(), "rejected by score:", rejected_score, "/", total_reads)
    
    # FINAL remaining BLAT-aligned queries:
    for query in query_details_dict:#query2gene_map:                        # contig/readID
        inner_dict = query_details_dict[query]
        if(inner_dict["match"]):
            mapped.add(query)
            db_match = inner_dict["gene"] 
            contig = inner_dict["is_contig"]
            bitscore = inner_dict["bitscore"]
            
            # RECORD alignment:
            if contig:                                      # If query is a contig, then
                for read in contig2read_map[query]:         #  take all reads making up that contig and
                    read_line = read + "<bitscore>" + str(bitscore)
                    if(db_match in gene2read_map):
                        gene2read_map[db_match].append(read_line)#  append their readIDs to aligned gene<->read dict,
                    else:
                        gene2read_map[db_match] = [read_line]
                        
            else:
                read_line = query + "<bitscore>" + str(bitscore)
                if(db_match in gene2read_map):
                    gene2read_map[db_match].append(read_line)
                else:
                    gene2read_map[db_match] = [read_line]

        
        #else:
        #    unmapped.add(query)
    return mapped, gene2read_map
    

def write_unmapped_seqs(unmapped_reads, reads_in, reads_out):
    #This is specifically made to export the reads as fastas, for the purpose of combining it with BLAT
    if(len(unmapped_reads) == 0):
        print(dt.today(), "no unmapped reads found.  skipping the write")
    else:
        read_seqs = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:])
        unmapped_seqs= []                                   # Initialize list of SeqRecords.
        for read in unmapped_reads:                         # Put corresponding SeqRecords for unmapped_reads
            unmapped_seqs.append(read_seqs[read])           #  into unmapped_seqs
        with open(reads_out,"w") as outfile:
            SeqIO.write(unmapped_seqs, outfile, "fasta")    #  and write it to file.

    # print no. aligned reads from current readtype set:
    #print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file) + '\n')
    #prev_mapping_count= len(mapped_reads)
    
def write_gene_map(DNA_DB, new_gene2read_file, gene2read_map, mapped_gene_file):
    # WRITE OUTPUT: rewrite gene<->read mapfile to include BLAT-aligned:
    # [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
    reads_count= 0
    genes= []
    with open(new_gene2read_file,"w") as out_map:               # Delete old gene2read_file and write a new one.
        for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                            #  (DNA db is needed to get the sequence.)
            if record.id in gene2read_map:                  #  If DNA db gene is one of the matched genes,
                genes.append(record)                        #  append the SeqRec to genes list (for next file), and
                out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                            #  write [aligned geneID, length, #reads, ...],
                for read in gene2read_map[record.id]:
                    out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                    reads_count+= 1
                else:
                    out_map.write("\n")                     #  and a new line character.

    # WRITE OUTPUT: BWA&BLAT-aligned geneIDs and seqs (.fna; fasta-format):
    with open(mapped_gene_file,"w") as outfile:
        SeqIO.write(genes, outfile, "fasta")    
        
       

def get_full_unmapped_reads(mapped_reads, fasta_keys):
    unmapped_reads = list()
    for item in fasta_keys:
        if (item not in mapped_reads):
            unmapped_reads.append(item)
    return unmapped_reads
    
    
def export_seqs(reads_in_dict, output_name):
    with open(output_name, "w") as out:
        for item in reads_in_dict:
            read_id = ">" + item + "\n"
            read_seq = reads_in_dict[item] + "\n"
            out.write(read_id)
            out.write(read_seq)
            

if __name__ == "__main__":

    DNA_DB              = sys.argv[1]   # INPUT: DNA db used for BLAT alignement
    contig_map_in       = sys.argv[2]   # INPUT: [contigID, #reads, readIDs ...]
    genes_out           = sys.argv[3]   # OUTPUT: BWA&BLAT-aligned geneIDs and aa seqs (.fna; fasta-format)
    gene_map_out        = sys.argv[4]   # OUTPUT: new gene2read_file instead of overwriting the old map
    
    reads_in            = sys.argv[5]   # INPUT: sequence in general
    blat_in             = sys.argv[6]
    reads_out           = sys.argv[7]


    input_safety = check_file_safety(reads_in) and check_file_safety(blat_in)
    
    if(input_safety):
        contig2read_map = import_contig_map(contig_map_in)
        blat_hits, reads_in_dict = get_blat_details(blat_in, reads_in)
        mapped_reads, gene2read_map = make_gene_map(blat_hits, contig2read_map) 
        unmapped_reads = get_full_unmapped_reads(mapped_reads, reads_in_dict)
        
        write_gene_map(DNA_DB, gene_map_out, gene2read_map, genes_out)
        write_unmapped_seqs(unmapped_reads, reads_in, reads_out)
        
    else:
        print(dt.today(), "input not safe.  either no reads, or empty BLATout.  copying")
        if(check_file_safety(reads_in)):
            copyfile(reads_in, reads_out)
    
