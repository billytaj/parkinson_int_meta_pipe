import os
import os.path
import sys
import pandas as pd


#This module takes in the Output report of the infernal tool, and the fastq it scanned.
#The goal is to bisect the fastq into 2 piles:  entries with IDs that were located by Infernal (rRNA)
#-> and entries that weren't (mRNA), and then export them

def extract_rRNA_ID(inf_file):
    #If there's a better way to do this, I'd like to see it.
    ID_list = set()
    inf_list = open(inf_file, mode='r')
    
    for item in inf_list:
        if(not item.startswith("#")):
            subitem_count = 0
            for subitem in item.split(' '):
                #our ID will always be the 3rd item, in a series of blanks.  
                #The blanks will always be there.  It's apparently not a tab
                if(len(subitem) > 0):
                    subitem_count += 1
                    print("subitem:", subitem, "length:", len(subitem))
                    if(subitem_count == 3):
                        #needs the @ at the start, for a match with the FASTQ's IDs
                        ID_list.add(str("@" + subitem))
                        break
                
        elif(len(item) == 2):
            break
    #pandas can't deal with sets, but we only need unique elements        
    return list(ID_list)

def filter_rRNA(rRNA_ID_list, fastq_sequence, mRNA_loc, rRNA_loc, file_name):
    #import the fastq as a DF
    fastq_df = pd.read_csv(fastq_sequence, header=None, names = [None])
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    #doing it inline so we don't create another DF
    #mRNA segment
    mRNA_export = os.path.join(mRNA_loc, file_name + "_mRNA.fastq")
    rRNA_export = os.path.join(rRNA_loc, file_name + "_rRNA.fastq")
    fastq_df[~fastq_df["ID"].isin(rRNA_ID_list)].to_csv(mRNA_export, sep = "\n", mode = "w+", header=False, index=False)
    fastq_df[fastq_df["ID"].isin(rRNA_ID_list)].to_csv(rRNA_export, sep="\n", mode = "w+", header=False, index=False)
    
    
if __name__ == "__main__":
    
    inf_file = sys.argv[1] #infernal
    fastq_sequence = sys.argv[2]
    mRNA_location = sys.argv[3]
    rRNA_location = sys.argv[4]
    segment_root_name = (inf_file.split('.')[0]).split("/")[-1]
    #segment root name expects the full path.  We're just reusing the name of the file, and changing the suffixes to our needs
    ID_list = extract_rRNA_ID(inf_file)
    filter_rRNA(ID_list, fastq_sequence, mRNA_location, rRNA_location, segment_root_name)
    