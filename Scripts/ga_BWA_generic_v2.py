#This script sifts through one SAM, and one reference to form a map of genes to constituent reads.
#Oct 22, 2020:
#new in v2: handle/discard multiple-hits to the same read
#couple of assumptions:
#-> contigs are now constructed of completely unique reads
#-> we do contigs now

import os
import os.path
import re
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from datetime import datetime as dt
import multiprocessing as mp
from shutil import copyfile

def get_match_score(cigar_segment):
    CIGAR = re.split("([MIDNSHPX=])", cigar_segment) # Split CIGAR string into list, placing
    CIGAR = CIGAR[:-1]                      #lop off the empty char artifact from the split
    position_count = 0                      #position counter, because the CIGAR string is split into alternating segments of <length><Label>, 
    length = 0
    matched = 0
    segment_length = 0
    for item in CIGAR:
        if((position_count %2) == 0):       #every even position (starting from 0) is going to be a length
            segment_length = int(item)
        elif((position_count %2) == 1):     #every odd position is going to be a label
            length += segment_length
            if(item == "M"):
                matched += segment_length
        position_count += 1
    if(length == 0):
        return 0
        
    match_score = 100 * (matched / length)
    
    return match_score

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


def import_contig2read(contig2read_file):
    # make initial dict of contigID<->readsID(s):
    contig2read_map = {}
    contig_reads = []                                    # list of just reads
    with open(contig2read_file,"r") as mapping:
        for line in mapping:
            if len(line)>5:                             # line starts with 'NODE_'
                entry= line.strip("\n").split("\t")     # break tab-separated into list
                contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
                contig_reads.extend(entry[2:])          # append all the reads
    
    return contig2read_map, contig_reads

def filter_common_contigs(contig2read_map, contig_reads):
    # make new dict only of contigs with unique reads:
    # (hard to tell w BWA what contigs match better, so for reads associated with multiple matched contigs, avoid choosing btw contigs for now.)
    contig_reads_count = Counter(contig_reads)           # dict of read<->no. of contigs
    contig2read_map_uniq = {}
    contig_unique_reads = []                             # DEBUG
    for contig in contig2read_map:
        for read in contig2read_map[contig]:            # If contig has
            if contig_reads_count[read]>1:              #  a read assoc. w multiple contigs
                break                                   #  then throw the contig away,
        else:
            contig2read_map_uniq[contig]= contig2read_map[contig]
                                                        #  else, store it in the unique dict.
            contig_unique_reads.extend(contig2read_map[contig]) # DEBUG

    return contig2read_map_uniq, contig_unique_reads

#####################################
# FUNCTION:
# add BWA-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:

# additional filtering steps:
# (1) only use contigs that contain unique reads (not shared w other contigs)
# (2) make sure BWA-aligned (.sam file has aligned & non-aligned data)
# (3) only use contigs/reads where >=90% length of the seq was matched

# CIGAR string describes how read aligns with the ref. Consists of >=1 components.
# Each component comprises an operator and no. bases which the op applies to.
#
# Operators:
# D	Deletion; the nucleotide is present in the reference but not in the read.
# H	Hard Clipping; the clipped nucleotides are not present in the read.
# I	Insertion; the nucleotide is present in the read  but not in the reference.
# M	Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
# N	Skipped region; a region of nucleotides is not present in the read.
# P	Padding; padded area in the read and not in the reference.
# S	Soft Clipping;  the clipped nucleotides are present in the read.
# =	Read Match; the nucleotide is present in the reference.
# X	Read Mismatch; the nucleotide is present in the reference.

def gene_map(sam, contig2read_map):#, mapped_reads, gene2read_map, contig2read_map):#, contig2read_map_uniq):                                      # Set of unmapped contig/readIDs=
                                                        #  gene_map(BWA .sam file)
    # tracking BWA-assigned & unassigned:
    query_details_dict = dict() #details about the match that we care about
    mapped = set()              #qualified mapped reads
    gene2read_map = dict()      #final gene->reads map
    unmapped = set()            #qualified unmapped reads
    
    len_chars= ["M","I","S","=","X"]                    # These particular CIGAR operations cause the
                                                        #  alignment to step along the query sequence.
                                                        # Sum of lengths of these ops=length of seq.

    # first, process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
            # "ON-THE-FLY" filter:
            inner_details_dict = dict() #stores the inner details: what's the gene association of this query(contig or read ID) and is it a contig+
            
            # extract & store data:
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  go to the next query (restart for).
            line_parts = line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and store:
            query = line_parts[0]                        #  queryID= contig/readID,
            db_match = line_parts[2]                     #  geneID, and a
            flag = bin(int(line_parts[1]))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                                                        #  where each bit is a flag for a specific descriptor.

            # is query BWA-aligned?:
            if flag[8]=="1":                            # If contig/read is NOT BWA-ALIGNED (9th digit=1),
                inner_details_dict["match"] = False
                query_details_dict[query] = inner_details_dict
                continue
                #unmapped.add(query)                     # add it to the unmapped set and
                #continue                                # skip to the next query.
            
            #if contig, mark it <we convert to reads down below>
            if query in contig2read_map:                
                contig = True                        
            else:
                contig = False                           
                
            # does query alignment meet threshold?:
            cigar_match_score = get_match_score(line_parts[5])
            if(cigar_match_score < 90):
                inner_details_dict["match"] = False
                query_details_dict[query] = inner_details_dict
                continue

            #reconcile multi-hit queries
            if(query in query_details_dict):
                if(query_details_dict[query]["match"]):
                    old_score = query_details_dict[query]["score"]
                    if(match_score > old_score):
                        #a better match came along.  if not, do nothing
                        query_details_dict[query]["score"] = match_score
                        query_details_dict[query]["gene"] = db_match
                        query_details_dict[query]["is_contig"] = contig
                else:
                    #previously unmatched
                    query_details_dict[query]["match"] = True
                    query_details_dict[query]["score"] = match_score
                    query_details_dict[query]["gene"] = db_match
                    query_details_dict[query]["is_contig"] = contig
            else:
                #new match
                inner_details_dict["match"] = True
                inner_details_dict["score"] = match_score
                inner_details_dict["gene"] = db_match
                inner_details_dict["is_contig"] = contig
                query_details_dict[query] = inner_details_dict
                
            
    
    for query in query_details_dict:
        inner_dict = query_details_dict[query]
        db_match = inner_dict["gene"]
        contig = inner_dict["is_contig"]    
        # RECORD alignments:
        if contig:                                      # If query is a contig, then
            gene2read_map[db_match].extend(contig2read_map[query])
                            
        else:                                           
            gene2read_map[db_match].append(query)       #  append its readID to aligned gene<->read dict,
        
        #sort the queries    
        if(inner_dict["match"]):
            mapped.add(query)
        else:
            unmapped.add(query)
            

    # return unmapped set:
    return unmapped, mapped, gene2read_map
    

def write_unmapped_reads(unmapped_reads, reads_in, output_file):
    if(len(unmapped_reads) == 0):
        print(dt.today(), "no unmapped reads found.  skipping")
    else:
        read_seqs = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:])
        # WRITE OUTPUT: non-BWA-aligned contig/readIDs:
        # and seqs (.fasta):
        unmapped_seqs = []                               # Inintialize list of SeqRecords.
        for read in unmapped_reads:                     # Put corresponding SeqRecords for unmapped_reads
            if(read in read_seqs):
                unmapped_seqs.append(read_seqs[read])       #  into unmapped_seqs
            else:
                print("ignoring:", read, "can't find in read_seqs")
        with open(output_file,"w") as out:
            SeqIO.write(unmapped_seqs, out, "fasta")    #  and write it to file.

        # print no. aligned reads from current readtype set:
        #print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(reads_in))
        #if x!=2: print ('')
        #prev_mapping_count= len(mapped_reads)


def write_gene_map(DNA_DB, gene2read_file, gene2read_map, mapped_gene_file):
    # WRITE OUTPUT: write gene<->read mapfile of BWA-aligned:
    # [BWA-aligned geneID, length, #reads, readIDs ...]
    reads_count = 0
    genes = []
    with open(gene2read_file,"w") as out_map:
        #btw, the "gene length" is literally the length of chars in the entry.  We can totally do away with seqIO
        for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                            #  (DNA db is needed to get the sequence.)
            if record.id in gene2read_map:                  #  If DNA db gene is one of the matched genes,
                genes.append(record)                        #  append the SeqRec to genes list (NOT REALLY USED), and
                out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                            #  write [aligned geneID, length, #reads, ...],
                for read in gene2read_map[record.id]:
                    out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                    reads_count+= 1
                else:
                    out_map.write("\n")                     #  and a new line character.
    
    #WRITE THE ANNOTATED GENES OUT TO A FILE.  FOR DOWNSTREAM USE
    with open(mapped_gene_file,"w") as outfile:
        SeqIO.write(genes, outfile, "fasta") 
    
if __name__ == "__main__":
    DNA_DB              = sys.argv[1]       # INPUT: DNA db used for BWA alignement
    contig2read_file    = sys.argv[2]       # INPUT: [contigID, #reads, readIDs ...]
    gene2read_out       = sys.argv[3]       # OUTPUT: [BWA-aligned geneID, length, #reads, readIDs ...]
    mapped_gene_file    = sys.argv[4]       # OUTPUT: genes mapped by BWA.
    
    reads_in            = sys.argv[5]   
    bwa_in              = sys.argv[6]
    reads_out           = sys.argv[7]
    
    input_safety = check_file_safety(reads_in) and check_file_safety(bwa_in)
    
    if(input_safety):
        contig2read_map, contig_reads = import_contig2read(contig2read_file)
        #contig2read_map_uniq, contig_unique_reads = filter_common_contigs(contig2read_map, contig_reads)
        # tracking BWA-assigned:
        gene2read_map = defaultdict(list)                    # dict of BWA-aligned geneID<->readID(s)
        mapped_reads = set()                                 # tracks BWA-assigned reads
        mapped_list = []
        prev_mapping_count = 0
        unmapped_reads, mapped_reads, gene2read_map = gene_map(bwa_in, mapped_reads, gene2read_map, contig2read_map, contig2read_map_uniq)
        
        write_gene_map(DNA_DB, gene2read_out, gene2read_map, mapped_gene_file)
        write_unmapped_reads(unmapped_reads, reads_in, reads_out)
        
    else:
        print(dt.today(), "input unsafe.  Either no reads, or BWA annotated nothing.  converting to fasta, then passing on")
        if(check_file_safety(reads_in)):
            if(reads_in.endswith(".fastq")):
                reads_to_convert = SeqIO.parse(reads_in, "fastq")
                SeqIO.write(reads_to_convert, reads_out, "fasta")
            else:
                copyfile(reads_in, reads_out)
            
    
    
    
    
    
    
    