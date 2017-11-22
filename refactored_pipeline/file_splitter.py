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
from Bio import SeqIO
import pandas as pd
import math as m

def split_fastq(file_name_in, dir_in, split_count = 4):
    #FASTQ has 4 lines per entry.
    file_base_name = file_name_in.split(".")[0]
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n")
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    #At this point, we've already got the number of reads.
    chunks = m.ceil(len(fastq_df) / split_count) 
    print("total df length:", len(fastq_df))
    print("chunk size:", chunks)
    if(chunks % 4) != 0:
        print("WARNING: split setting will not yield even divisions")
        while((chunks * split_count) < int(len(fastq_df))):
            chunks += 1
        print("new compensatory chunk size:", chunks)
    #sys.exit()
    for i in range(0, split_count):
        print("working on segment :", i+1, "of", split_count)
        if not os.path.exists(dir_in):
            os.makedirs(dir_in)
        new_file_name = dir_in + file_base_name + "_"+str(i) + ".fastq"
        start_index = int(i * chunks)
        end_index = int(((i+1) * chunks)-1)
        print("[", i, "]: start:", start_index, "| end:", end_index)
        fastq_df.iloc[start_index:end_index, :].to_csv(new_file_name, chunksize = chunks, index=False, sep='\n', header=False)

    
if __name__ == "__main__":
    input_file = sys.argv[1]
    
    input_extension = input_file.split(".")[1]
    operating_dir = "C:\\parkinson_int_meta_pipe\\refactored_pipeline\\"
    os.chdir(operating_dir)
    if(input_extension == "fastq"):
        #print("This is a fastq")    
        split_fastq(input_file, operating_dir + "\\new_dump\\", 3)
    else:
        print("This is not a fasta.  Not splitting the file")
    sys.exit()