#!/usr/bin/env python
# Oct 16, 2017
# -------------------------------------------------
# this module looks to be responsible for splitting up the fasta or fastq.

# Nov 21, 2017
# -------------------------------------- 
# upon further inspection, this thing really does just split up the fastq, and fasta.  but with unnecessary logic
# According to BJ:  splitting fasta isn't necessary anymore.  It was a leftover piece of code.  
# This thing will only split Fastq until we see a need for something else

import sys
import os
import os.path
import shutil
import subprocess
import multiprocessing
#from Bio import SeqIO
import pandas as pd
import math as m

def split_fastq(file_name_in, file_name_out, split_count = 4):
    #FASTQ has 4 lines per entry.
    file_base_name = file_name_in.split(".")[0]
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n")
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    #At this point, we've already got the number of reads.
    chunks = m.ceil(len(fastq_df) / split_count) #how many sequences each split file will have
    print("total df length:", len(fastq_df))
    print("chunk size:", chunks)
    if(chunks < 1):
        print("split count too large. not enough info to split")
        chunks = 1
        print("new split count:", split_count)
        
    elif(chunks % 4) != 0:
        print("WARNING: split setting will not yield even divisions")
        while((chunks * split_count) < int(len(fastq_df))):
            chunks += 1
        print("new compensatory chunk size:", chunks)
        
    for i in range(0, split_count):
        print("working on segment :", i+1, "of", split_count)
        #fancy naming
        new_file_name = file_name_out + "_" + str(i) + ".fastq"
        if(split_count == 1):
            new_file_name = file_name_out + ".fastq"
        #split file by selective selection, and writing
        start_index = int(i * chunks)
        end_index = int(((i+1) * chunks))
        #if(chunks == 1):
        #    end_index += 1 #override on splits that only have 1 
        if not(fastq_df.iloc[start_index:end_index, :].empty):
            fastq_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize = chunks, mode = "w+", index=False, sep='\n', header=False)
        else:
            print("empty frame detected.  no sense in running the rest")
            break

    
if __name__ == "__main__":
    if(len(sys.argv) == 4):
        input_file = sys.argv[1]
        output_name = sys.argv[2]
        split_count = int(sys.argv[3])
        if(split_count < 2):
            print("not splitting file, just moving it")
            split_count = 1
            
        input_extension = input_file.split(".")[1]
        if(input_extension == "fastq"):
            #print("This is a fastq")    
            split_fastq(input_file, output_name, split_count)
        else:
            print("This is not a fastq.  Not splitting the file")
    else:
        print("only: ", len(sys.argv), "number of args.  not running")
        print("wrong number of args to file splitter.  not running")
    sys.exit()