import os
import sys
import pandas

# This piece of code was originally based off of the map_read_gene_BWA.
# Its purpose is filter out what's already been scanned inside the various sequence sections, from the contig assembly.
# So, the 
# This code has 3 parts:
# 1) Import the contig manifest
# 2) import the non-redundant microbial cds database from NIH (comes in a FASTA format)
# 

def pass_quality(cigar):
    # This returns the quality of the match, as a percentage.  
    length= 0
    matched= 0
    len_chars= ["M","I","S","=","X"]  
    CIGAR= re.split("([MIDNSHPX=])", cigar) # Split CIGAR string into list, placing
                                                    #  all chars w/in [...] into own field
                                                    #  (e.g., 9S41M50S->['9','S','41','M','50','S','']).
    for index in range(len(CIGAR))[:-1]:            # Loop CIGAR elements (last element=''),
        if CIGAR[index+1] in len_chars:             # Use CIGAR operations that step along the query seq,
            length+= int(CIGAR[index])              #  to determine length of query.
        if CIGAR[index+1]=="M":                     # Use CIGAR match operation to
            matched+= int(CIGAR[index])
    if(length > 0):        
        return matched / length
    else:
        return 0

def make_gene_map():
    pair_1_sam_df = pd.read_csv(sam_file, error_bad_lines=False, header=None, sep="\t")                 # read it
    pair_1_sam_df.iloc[:, 1] = pair_1_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])   # rewrite the CIGAR flag
    selected_mapped = pair_1_sam_df.loc[pair_1_sam_df.iloc[:, 5] < "1"].iloc[:, 0:6]                   # select what we need from the SAM.  unmapped reads only
    selected_mapped.columns = ["reads", "flag", "contig"]                                               # name the columns for easy manipulation 
    selected_mapped.index = selected_mapped.groupby('contig').cumcount()                                # this looks like it makes the indices, by finding the longest row
    t2 = selected_mapped.pivot(values='reads', columns='contig')                                        # pivot it so the contig is the columns
    new_t2 = t2.T   
    
    return new_t2
    
if __name__ == "__main__":
    contig_manifest = sys.argv[1] # the contig2read file, which is a map that contains contig, the number of reads, and the inidividual reads that constitute the contig
    