#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - added .strip("\n") to extracting from the gene2read_file
# - fixed sort function to return float
# - reversed the sort direction to get largest scores at top
# - pulled gene_map() out of the readtype loop
# - rewrote gene_map()
# - removed checks to see if contig/read was already assigned to same gene (shouldn't happen)
# - removed check to see if query=read was already mapped by BLAT (shouldn't happen)
# - changed align_len to float(align_len) (to make sure it's a number)
# - paired-end BLAT outputs are combined and analyzed at the same time

# NOTE:
# - Filenames for unmerged paired-end reads must be specified last.
# - Any read aligned by BWA will not be input into BLAT. Even BWA-aligned contigs
#   have unique reads (not found in other contigs).
# - BLATpp only takes the read<->gene match with the top score. For duplicate
#   top scores, it only takes the one that shows up first in the initial sort.
# - The paired-end BLAT outputs are combined and analyzed at the same time, in order to
#   rank the matches between different ends of the same read, and not double count that read.
# - Sometimes a read is aligned via multiple contigs and annoteted only according
#   to the best matched contig only. Theoretically, this breaks up the other contigs.

import os
import os.path
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from datetime import datetime as dt
import multiprocessing as mp

def import_contig_map(contig2read_file):
    # make dict of contigID<->readsID(s):
    contig2read_map= {}
    with open(contig2read_file,"r") as mapping:
        for line in mapping:
            if len(line)>5:                             # line starts with 'NODE_'
                entry= line.strip("\n").split("\t")     # break tab-separated into list
                contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
    return contig2read_map


def import_gene_map(gene2read_file):
    # make dict of BWA-aligned geneID<->readID(s):
    BWAreads = []                                        # DEBUG
    gene2read_map = defaultdict(list)                    # Dict of BWA&BLAT-aligned geneID<->readID(s)
    mapped_reads = set()
    with open(gene2read_file,"r") as mapping:           #  initialized w BWA-alignments.
        for line in mapping:
            if len(line)>5:                             # line at least 5 characeters?
                entry = line.strip("\n").split("\t")
                gene2read_map[entry[0]] = entry[3:]      # key=geneID, value=list of readID(s)
                                                        # (using [:] syntax ensures a list, even if one ele)
                mapped_reads.update(entry[3:])
                BWAreads.extend(entry[3:])              # DEBUG
    return gene2read_map, mapped_reads, BWAreads
    

#####################################
# FUNCTIONS:

# read .blatout file and acquire read length:
def get_blat_details(blat_in, reads_in):                          # List of lists [blatout info, read length]=
                                                        #  read_aligned(.blatout file, dict contig/readID<->SeqRecord)
    seqrec = SeqIO.index(reads_in, os.path.splitext(reads_in)[1][1:])
    # get info from .blatout file:
    print ('Reading ' + str(os.path.basename(blat_in)) + '.')
    with open(blat_in,"r") as tabfile:
        hits= []                                        # List of lists containing .blatout fields.
        for line in tabfile:                            # In the .blatout file:
            if len(line)>=2:                            # If length of line >= 2,
                info_list= line.strip("\n").split("\t") #  make a list from .blatout tab-delimited fields,
                query= info_list[0]                     #  get the queryID= contig/readID,
                info_list.append(len(seqrec[query].seq))#  add query length (int) to end of list,
                hits.append(info_list)                  #  and append the list to the hits list.

    # return info:
    print(dt.today(), "hits in blatout:", len(hits))
    return hits

# sort by score:
# (12th field of the .blatout file)
def sortbyscore(line):
    return float(line[11])

# add BLAT-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:
def gene_map(hits, mapped_reads, gene2read_map, contig2read_map):                         # fail-mapped contig/readIDs=
                                            #  gene_map(list of list of blatout fields)
    print(dt.today(), "number of keys in gene map [prior]:", len(gene2read_map.keys())) 
    print(dt.today(), "number of mapped reads [prior]:", len(mapped_reads))
    # sort alignment list by high score:
    sorted_hits= sorted(hits, key=sortbyscore, reverse=True)
    del hits

    # BLAT threshold:
    identity_cutoff= 85
    length_cutoff= 0.65
    score_cutoff= 60

    # tracking BLAT-assigned & unassigned:
    #query2gene_map= defaultdict(set)        # Dict of BLAT-aligned contig/readID<->geneID(s).
    #queryIScontig= {}                       # Dict of BLAT-aligned contig/readID<->contig? (True/False).
    unmapped = set()                         # Set of unmapped contig/readIDs.
    query_details_dict = dict()

    # loop through sorted BLAT-aligned reads/contigs:
    for line in sorted_hits:
        inner_details_dict = dict()
        # "ON-THE-FLY" filter:
        
        # extract & store data:
        query = line[0]                      # queryID= contig/readID
        db_match = line[1]                   # geneID
        seq_identity = float(line[2])        # sequence identity
        align_len = int(line[3])             # alignment length
        score = float(line[11])              # score
        seq_len = int(line[12])              # query sequence length
        
        # is this alignment the highest-score match for that query?
        if query in query_details_dict: #query2gene_map:         # If alignment previously found for this contig/read,
            continue                        # skip to next qurey as this alignment will have a lower score,
                                            # and don't add to unmapped set.
        # is query a contig?:
        if query in contig2read_map:        # If query is found in contig list,
            contig= True                    #  then mark as contig,
        else:                               #  if not,
            contig= False                   #  mark as not a contig.
                                        
        # does query alignment meet threshold?:
        # (identity, length, score):
        if seq_identity<=identity_cutoff or align_len<=seq_len*length_cutoff or score<=score_cutoff:
            unmapped.add(query)             # if threshold is failed, add to unmapped set and
            continue                        # skip to the next query.
        
        # store info for queries that remain:
        inner_details_dict["is_contig"] = contig
        inner_details_dict["gene"] = db_match
        query_details_dict[query] = inner_details_dict
        #queryIScontig[query] = contig        # Store contig (T/F) info.
        #query2gene_map[query].add(db_match) # Collect all aligned genes for contig/read;
                                            #  query2gene_map[query] is a set, although there should be
                                            #  no more than one gene alignement getting through filter.
        
    # EXPAND HERE to deal with multiple high-score genes for a read.
    
    # DEBUG:
    contigread_inmapped = 0
    contigread_inmapped_ingene = 0
    contigread_ingene = 0
    read_inmapped = 0
    read_inmapped_ingene = 0
    read_ingene = 0
    
    # FINAL remaining BLAT-aligned queries:
    for query in query_details_dict:#query2gene_map:                        # contig/readID
        inner_dict = query_details_dict[query]
        db_match = inner_dict["gene"] #list(query2gene_map[query])[0]        # geneID (pull out of 1-element set)
        contig = inner_dict["is_contig"] #queryIScontig[query]                    # contig?
        
        # RECORD alignment:
        if contig:                                      # If query is a contig, then
            for read in contig2read_map[query]:         #  take all reads making up that contig and
            
                # DEBUG:
                if read in mapped_reads:                # (Track how many contig reads have already been
                    contigread_inmapped+= 1             #  mapped by BLAT---i.e., through other contigs---
                    if read in gene2read_map[db_match]: #  and how many of those were already assigned to this
                        contigread_inmapped_ingene+= 1  #  particular gene---i.e., contigs aligned to same gene.)
                elif read in gene2read_map[db_match]:   # (Check to see if any of the contig reads are already
                    contigread_ingene+= 1               #  assigned to this particular gene, but not by BLAT.)
                
                if read not in mapped_reads:            #  if not already assigned by BLAT to a diff gene***, then
                    gene2read_map[db_match].append(read)#  append their readIDs to aligned gene<->read dict,
                    mapped_reads.add(read)              #  and mark them as assigned by BLAT.
        elif not contig:                                # If query is a read, then
        
            # DEBUG:
            if query in mapped_reads:                   # (Check to see if the read has already been mapped
                read_inmapped+= 1                       #  by BLAT---it shouldn't be---and
                if query in gene2read_map[db_match]:    #  and to the same gene, for that matter.
                    read_inmapped_ingene+= 1
            elif query in gene2read_map[db_match]:      # (Check to see if the read is already assigned to this
                read_ingene+= 1                         #  particular gene, but not by BLAT---it shouldn't be.)
            
            if query not in mapped_reads:
                gene2read_map[db_match].append(query)       #  append its readID to aligned gene<->read dict,
                mapped_reads.add(query)                     #  and mark it as assigned by BLAT.

        # *** This deals with reads that show up in multiple contigs.
        # Just use the read alignment from the contig that had the best alignment score.
        # This could result in "broken" contigs...

    # DEBUG (for this datatype):
    print ('no. contig reads already mapped by BLAT = ' + str(contigread_inmapped))
    print ('no. contig reads mapped by BLAT to same gene = ' + str(contigread_inmapped_ingene))
    print ('no. contig reads mapped by NOT BLAT to same gene = ' + str(contigread_ingene) + ' (should be 0)')
    print ('no. reads already mapped by BLAT = ' + str(read_inmapped) + ' (should be 0)')
    print ('no. reads already mapped by BLAT to same gene = ' + str(read_inmapped_ingene) + ' (should be 0)')
    print ('no. reads already mapped by NOT BLAT to same gene = ' + str(read_ingene) + ' (should be 0)')

    # Remove contigs/reads previously added to the unmapped set but later found to have a mapping:
    # This prevents re-annotation by a later program.
    # Such queries end up in the unmapped set when they BLAT-aligned to multiple genes, where one
    # alignment is recorded, while the other alignments fail the "on-the-fly" alignment-threshold filter.
    print ('umapped no. (before double-checking mapped set) = ' + str(len(unmapped)))
    for query in query_details_dict: #query2gene_map:                        # Take all contigs/reads to be mapped and
        try:                                            #  if they exist in the unmapped set, then
            unmapped.remove(query)                      #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))

    # return unmapped set:
    print(dt.today(), "number of keys in gene map [post]:", len(gene2read_map.keys())) 
    print(dt.today(), "number of mapped reads [post]:", len(mapped_reads))
    
    
    return unmapped    
def write_unmapped_seqs(unmapped_reads, reads_in, reads_out):
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
    # (this wasn't done in BWA post-processing)
    with open(mapped_gene_file,"w") as outfile:
        SeqIO.write(genes, outfile, "fasta")    

if __name__ == "__main__":

    DNA_DB= sys.argv[1]                 # INPUT: DNA db used for BLAT alignement
    contig2read_file= sys.argv[2]       # INPUT: [contigID, #reads, readIDs ...]
    gene2read_file= sys.argv[3]         # INPUT: [BWA-aligned geneID, length, #reads, readIDs ...]
                                        # ->OUTPUT: [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
    mapped_gene_file= sys.argv[4]       # OUTPUT: BWA&BLAT-aligned geneIDs and aa seqs (.fna; fasta-format)
    new_gene2read_file = sys.argv[5]    # OUTPUT: new gene2read_file instead of overwriting the old map
    new_gene_file = sys.argv[6]         # OUTPUT new gene.fna 

    operating_mode = "single"
    contigs_reads_in = sys.argv[7]
    contigs_blat_in = sys.argv[8]
    contigs_reads_out = sys.argv[9]

    singletons_reads_in = sys.argv[10]
    singletons_blat_in = sys.argv[11]
    singletons_reads_out = sys.argv[12]

    if(len(sys.argv) == 19):
        operating_mode = "paired"
        pair_1_reads_in = sys.argv[13]
        pair_1_blat_in = sys.argv[14]
        pair_1_reads_out = sys.argv[15]
        
        pair_2_reads_in = sys.argv[16]
        pair_2_blat_in = sys.argv[17]
        pair_2_reads_out = sys.argv[18]
        print(dt.today(), "OPERATING MODE:", operating_mode)
    elif(len(sys.argv) == 13):
        print(dt.today(), "OPERATING MODE:", operating_mode)
    else:
        print(dt.today(), "wrong number of args.  something wrong in metapro_commands. exiting")
        sys.exit()

    contig2read_map = import_contig_map(contig2read_file)
    # tracking BLAT-assigned:
    prev_mapping_count= 0
    gene2read_map, mapped_reads, BWAreads = import_gene_map(gene2read_file)

    # DEBUG:
    if len(set(BWAreads))==len(BWAreads):
        print ('BWA-aligned reads are all unique.\n')
    else:
        print ('There are repeating BWA-aligned reads:')
        print ('no. unique reads= ' + str(len(set(BWAreads))))
        print ('no. total reads= ' + str(len(BWAreads)) + '\n')
        print(dt.today(), "THIS IS A PROBLEM. shutting down")
        sys.exit()

    

    contigs_blat_hits = get_blat_details(contigs_blat_in, contigs_reads_in)
    contigs_unmapped_reads = gene_map(contigs_blat_hits, mapped_reads, gene2read_map, contig2read_map)

    singletons_blat_hits = get_blat_details(singletons_blat_in, singletons_reads_in)
    singletons_unmapped_reads = gene_map(singletons_blat_hits, mapped_reads, gene2read_map, contig2read_map)

    if(operating_mode == "paired"):
        pair_1_blat_hits = get_blat_details(pair_1_blat_in, pair_1_reads_in)
        pair_1_unmapped_reads = gene_map(pair_1_blat_hits, mapped_reads, gene2read_map, contig2read_map)
        
        #pair_2_blat_hits = get_blat_details(pair_2_blat_in, pair_2_reads_in)
        #pair_2_unmapped_reads = gene_map(pair_2_blat_hits, mapped_reads, gene2read_map, contig2read_map)

    process_store = []  

    gene_map_write_process = mp.Process(target = write_gene_map, args = (DNA_DB, new_gene2read_file, gene2read_map, mapped_gene_file))
    gene_map_write_process.start()
    print(dt.today(), "GA BLAT gene map export launched")
    process_store.append(gene_map_write_process)
    
    contig_write_unmapped_process = mp.Process(target = write_unmapped_seqs, args = (contigs_unmapped_reads, contigs_reads_in, contigs_reads_out))
    contig_write_unmapped_process.start()
    print(dt.today(), "GA BLAT unmapped contigs export launched")
    process_store.append(contig_write_unmapped_process)

    singletons_write_unmapped_process = mp.Process(target = write_unmapped_seqs, args = (singletons_unmapped_reads, singletons_reads_in, singletons_reads_out))
    singletons_write_unmapped_process.start()
    print(dt.today(), "GA BLAT unmapped singletons export launched")
    process_store.append(singletons_write_unmapped_process)

    if(operating_mode == "paired"):
        pair_1_write_unmapped_process = mp.Process(target = write_unmapped_seqs, args = (pair_1_unmapped_reads, pair_1_reads_in, pair_1_reads_out))
        pair_1_write_unmapped_process.start()
        print(dt.today(), "GA BLAT unmapped pair 1 export launched")
        process_store.append(pair_1_write_unmapped_process)
        
        pair_2_write_unmapped_process = mp.Process(target = write_unmapped_seqs, args = (pair_1_unmapped_reads, pair_2_reads_in, pair_2_reads_out))
        pair_2_write_unmapped_process.start()
        print(dt.today(), "GA BLAT unmapped pair 2 export launched")
        process_store.append(pair_2_write_unmapped_process)

    for item in process_store:
        item.join()
    process_store[:] = []
    print(dt.today(), "GA BLAT finished")

