import os
import sys
import time
from datetime import datetime as dt

#this code will merge all FASTAs with the same source 
#it's meant to condense all of the leftover FASTA reads from running them in different segments of chocophlan via BWA
#notes: overlaps in IDs are gauranteed to be the same data.  So we don't even need to check if a read ID has been used before

if __name__ == "__main__":
    in_dir = sys.argv[1] #eg: GA_BWA/
    segment_name = sys.argv[2]  #eg: pair_1_00, contigs_1, singletons_9111
    out_dir = sys.argv[3]
    out_path = os.path.join(out_dir, segment_name + ".fasta")
    dir_list = os.listdir(in_dir)
    

    #merge it all into a dict.  Let the duplicate IDs overwrite. it's the same data
    read_dict = dict()
    for item in dir_list:
        if(item.startswith(segment_name)):
            file_to_open = os.path.join(in_dir, item)
            with open(file_to_open, "r") as in_fasta:
                line_is_ID = True
                ID = ""
                seq = ""
                for line in in_fasta:
                    if(line_is_ID):
                        ID = line.strip("\n")
                        line_is_ID = False
                    else:
                        seq = line.strip("\n")
                        line_is_ID = True
                        read_dict[ID] = seq

    #export it                    
    with open(out_path, "w") as out_fasta:
        for item in read_dict:
            out_line = item + "\n"
            out_fasta.write(out_line)
            out_line = read_dict[item] + "\n"
            out_fasta.write(out_line)
    