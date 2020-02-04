#This script will take in the BWA SAM file from contigs, the gene prediction log from 
#GeneMark, and assign the reads to the segments, based on the alignment position within the contig
#This importer no longer needs the SAM to be trimmed.  
import os
import sys
import pandas as pd

def import_bwa_sam(bwa_whiole_contig_sam):
    with open(bwa_whole_contig_sam, "r") as bwa_sam:
        for line in bwa_sam:
            if(line.startswith("@PG")):
                continue
            elif(line.startswith("@SQ")):
                continue
            elif(len(line) <= 1):
                continue
            else:
                sam_align = line.strip("\n").split("\t")
                read_ID = sam_align[0]
                flag = bin(int(sam_align[1]))[2:].zfill(11)
                if(flag[8] == "1"):
                    #read is unmapped
                    print(read_ID, "is unmapped")
                    #continue
                else:
                    print(read_ID, "found")
                    

if __name__ == "__main__":
    bwa_sam = sys.argv[1]
    import_bwa_sam(bwa_sam)