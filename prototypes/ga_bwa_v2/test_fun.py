import os
import os.path
import re
import sys
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
import pandas as pd


def gene_map(sam):                                      # Set of unmapped contig/readIDs=
                                                        #  gene_map(BWA .sam file)
    # tracking BWA-assigned & unassigned:
    query2gene_map= defaultdict(set)                    # Dict of BWA-aligned contig/readID<->geneID(s).
    queryIScontig= {}                                   # Dict of BWA-aligned contig/readID<->contig? (True/False).
    unmapped= set()                                     # Set of unmapped contig/readIDs.
    mapped_reads = set()
    mapped_list = []

    len_chars= ["M","I","S","=","X"]                    # These particular CIGAR operations cause the
                                                        #  alignment to step along the query sequence.
                                                        # Sum of lengths of these ops=length of seq.

    # first, process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
        
            # "ON-THE-FLY" filter:
        
            # extract & store data:
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  go to the next query (restart for).
            line_parts= line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and store:
            query= line_parts[0]                        #  queryID= contig/readID,
            db_match= line_parts[2]                     #  geneID, and a
            flag= bin(int(line_parts[1]))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                                                        #  where each bit is a flag for a specific descriptor.

            # is query BWA-aligned?:
            if flag[8]=="1":                            # If contig/read is NOT BWA-ALIGNED (9th digit=1),
                unmapped.add(query)                     # add it to the unmapped set and
                continue                                # skip to the next query.
            
            # is query a contig made of unique reads?:
            if query in contig2read_map:                # If query is a contig (searches through keys)
                if query in contig2read_map_uniq:       #  and it's made of contig-unique reads,
                    contig= True                        #  then mark as contig and move on.
                else:                                   # Otherwise, contig contains non-unique reads,
                    unmapped.add(query)                 #  therefore add it to the unmapped set and
                    continue                            #  skip to the next query.
            else:
                contig= False                           # If query isn't a contig, just move on.
            
            # does query alignment meet threshold?:
            length= 0
            matched= 0
            CIGAR= re.split("([MIDNSHPX=])", line_parts[5]) # Split CIGAR string into list, placing
                                                            #  all chars w/in [...] into own field
                                                            #  (e.g., 9S41M50S->['9','S','41','M','50','S','']).
            for index in range(len(CIGAR))[:-1]:            # Loop CIGAR elements (last element=''),
                if CIGAR[index+1] in len_chars:             # Use CIGAR operations that step along the query seq,
                    length+= int(CIGAR[index])              #  to determine length of query.
                if CIGAR[index+1]=="M":                     # Use CIGAR match operation to
                    matched+= int(CIGAR[index])             #  determine no. nuclotides matched.
            if matched<length*0.9:                          # If alignment is <90% matched:
                unmapped.add(query)                         # add it to the unmapped set and
                continue                                    # skip to the next query.

            # store info for queries that remain:
            query2gene_map[query].add(db_match)             # Collect all aligned genes for contig/read.
            queryIScontig[query]= contig                    # Store contig (T/F) info.

    print ('Reading ' + str(os.path.basename(sam)) + '.')
    print ('no. queries (that meet initial threshold)= ' + str(len(query2gene_map)))
    
    # remove read queries that are part of contigs:
    # THOUGH THIS SHOULDN'T HAPPEN (and it doesn't).
    for query in query2gene_map:                        # Check all queries
        if queryIScontig[query]:                        #  (that are reads) to see
            if query in contig_reads:                   #  if they belong to any contigs, and if so
                del query2gene_map[query]               #  delete the read from the query list.
                del queryIScontig[query]                # Don't add it to the unmapped set, as it shouldn't be
                                                        #  annotated as a stand-alone read.

    print ('no. (after removal of read queries in contigs)= ' + str(len(query2gene_map)))

    # remove queries aligning to multiple genes:
    for query, genes in list(query2gene_map.items()):   # Check all queries to see
        if len(genes)>1:                                #  if they are aligend to multiple genes, and if so
            unmapped.add(query)                         #  add it to the unmapped set, then
            del query2gene_map[query]                   #  delete it from the query list.
            del queryIScontig[query]

    print ('no. (after removal of queries aligning to multiple genes)= ' + str(len(query2gene_map)))

    # FINAL remaining BWA-aligned queries:
    for query in query2gene_map:                        # contig/readID
    
        db_match= list(query2gene_map[query])[0]        # geneID (pull out of 1-element set)
        contig= queryIScontig[query]                    # contig?
    
        # RECORD alignments:
        if contig:                                      # If query is a contig, then
            gene2read_map[db_match].extend(contig2read_map_uniq[query])
                                                        #  add the readIDs of all the reads making up
                                                        #  that contig to the aligned gene<->read dict,
            for read in contig2read_map_uniq[query]:    #  and mark all those reads
                mapped_reads.add(read)                  #  as assigned by BWA.
                mapped_list.append(read)                # DEBUG (not a set so will double w paired reads)
        else:                                           # If query is a read,
            gene2read_map[db_match].append(query)       #  append its readID to aligned gene<->read dict,
            mapped_reads.add(query)                     #  and mark it as assigned by BWA.
            mapped_list.append(query)                   # DEBUG (not a set so will double w paired reads)

    # Remove contigs/reads previously added to the unmapped set but later found to have a mapping:
    # This prevents re-annotation by a later program.
    # Such queries end up in the unmapped set when they BWA-aligned to multiple genes, where one
    # alignment is recorded, while the other alignments fail the "on-the-fly" filter; or
    # when one end of an unmerged paired read aligns and the other end doesnt.
    print ('umapped no. (before double-checking mapped set)= ' + str(len(unmapped)))
    for query in query2gene_map:                        # Take all contigs/reads to be mapped and
        try:                                            #  if they exist in the unmapped set, then
            unmapped.remove(query)                      #  remove them from the unmapped set.
        except:
            pass
    print ('umapped no. (after double-checking mapped set)= ' + str(len(unmapped)))


    # DEBUG (check so far):
    if len(set(mapped_list))==len(mapped_list):
        print ('So far, BWA-aligned reads are all unique.')
    else:
        print ('Repeating BWA-aligned reads thus far:')
        print ('no. unique reads= ' + str(len(set(mapped_list))))
        print ('no. total reads= ' + str(len(mapped_list)))
        print ('no reads in the set= ' + str(len(mapped_reads)))

    # return unmapped set:
    return (unmapped, mapped_reads, mapped_list)

#####################################

def import_fasta(file_name_in):
    #this imports the fasta into a pandas dataframe
    fasta_df = pd.read_csv(file_name_in, error_bad_lines=False, header=None, sep="\n")  # import the fasta
    fasta_df.columns = ["row"]
    #There's apparently a possibility for NaNs to be introduced in the raw fasta.  We have to strip it before we process (from DIAMOND proteins.faa)
    fasta_df.dropna(inplace=True)
    new_df = pd.DataFrame(fasta_df.loc[fasta_df.row.str.contains('>')])  # grab all the IDs
    new_df.columns = ["names"]
    new_data_df = fasta_df.loc[~fasta_df.row.str.contains('>')]  # grab the data
    new_data_df.columns = ["data"]
    fasta_df = new_df.join(new_data_df, how='outer')  # join them into a 2-col DF
    fasta_df["names"] = fasta_df.fillna(method='ffill')  # fill in the blank spaces in the name section
    fasta_df["names"] = fasta_df["names"].apply(lambda x: x.split(" ")[0])
    fasta_df.dropna(inplace=True)  # remove all rows with no sequences
    fasta_df.index = fasta_df.groupby('names').cumcount()  # index it for transform
    temp_columns = fasta_df.index  # save the index names for later
    fasta_df = fasta_df.pivot(values='data', columns='names')  # pivot
    fasta_df = fasta_df.T  # transpose
    fasta_df["sequence"] = fasta_df[fasta_df.columns[:]].apply(lambda x: "".join(x.dropna()), axis=1)  # consolidate all cols into a single sequence
    fasta_df.drop(temp_columns, axis=1, inplace=True)
    #fasta_df["names"] = fasta_df.index
    return fasta_df


def filter_consumed_reads(read_seqs, BWA_sam_file, output_file, is_unmerged_flag):
    #if (not os.path.exists(read_seqs)):
    #    return None
    if (not os.path.exists(BWA_sam_file)):
        return None
    #if numsets not in [2,4]:
    #    sys.exit('Incorrect number of readtype sets.')

    # process BWA output:
    # readtype sets: contigs, merged, unmerged1, unmerged2
    #for x in range(int((len(sys.argv)-4)/3)):
    #read_file= sys.argv[3*x+4]      # INPUT: all contig/readIDs and seqs (.fasta)
    #read_seqs= SeqIO.index(read_file, os.path.splitext(read_file)[1][1:]) #This line just grabs all the IDs from the input fasta
                                    # dict of all read SeqRecords: key=contig/readID
                                    #  (second argument specifies filetype, e.g., "fasta")
    #BWA_sam_file= sys.argv[3*x+5]   # INPUT: BWA-aligned&unaligned contig/readIDs (.sam)
    #output_file= sys.argv[3*x+6]    # OUTPUT: non-BWA-aligned contig/readIDs and seqs (.fasta)
    #read_seqs = 
    count = 0
    for item in read_seqs:
        print(item)
        count += 1
        if(count > 10):
            break

    # extraction of "non-BWA-aligned" and "BWA-aligned":
    unmapped_reads = []
    if (is_unmerged_flag):                          # Only do once for unmerged paired end (same .sam file).
        unmapped_reads, mapped_reads, mapped_list = gene_map(BWA_sam_file)      # Store BWA-aligned contigs/reads in gene2read_map
                                                    #  (aligned geneID<->readID(s) dict),
                                                    #  and return a set of unmapped readIDs
    # WRITE OUTPUT: non-BWA-aligned contig/readIDs:
    # and seqs (.fasta):
    unmapped_seqs= []                               # Inintialize list of SeqRecords.
    for read in unmapped_reads:                     # Put corresponding SeqRecords for unmapped_reads
        if(read in read_seqs):
            unmapped_seqs.append(read_seqs[read])       #  into unmapped_seqs
        else:
            print("ignoring:", read, "can't find in read_seqs")
    with open(output_file,"w") as out:
        SeqIO.write(unmapped_seqs, out, "fasta")    #  and write it to file.

    # print no. aligned reads from current readtype set:
    #print (str(len(mapped_reads)-prev_mapping_count) + ' additional reads were mapped from ' + os.path.basename(read_file))
    #prev_mapping_count= len(mapped_reads)

if __name__ == "__main__":
    contig_file = sys.argv[1]
    sam_file = sys.argv[2]
    output_file = sys.argv[3]
    #umapped, mapped_reads, mapped_list = gene_map(sam_file)
    contig_df = import_fasta(contig_file)
    filter_consumed_reads(contig_df, sam_file, output_file, True)
