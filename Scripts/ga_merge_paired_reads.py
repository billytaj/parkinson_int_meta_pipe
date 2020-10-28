import os
import sys

import pandas as pd
#Oct 27, 2020
#This code merges the paired reads together.  It concatenates p1, and the reverse-complement of p2 to form a super-read.


def import_fastq(file_name_in):
    fastq_df = pd.read_csv(file_name_in, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "sequences", "junk", "quality"]
    fastq_df["ID"] = fastq_df["ID"].apply(lambda x: x.strip("@"))
    return fastq_df

def reverse_complement_seq(x):
    final_string = ""
    rev = x[::-1]
    for item in rev:
        if(item == "T"):
            final_string += "A"
        elif(item =="G"):
            final_string += "C"
        elif(item == "A"):
            final_string += "T"
        elif(item == "C"):
            final_string += "G"
    return final_string

def reverse_quality(x):
    return x[::-1]

if __name__ == "__main__":
    pair_1_file = sys.argv[1]
    pair_2_file = sys.argv[2]
    export_file = sys.argv[3]
    pair_1_df = import_fastq(pair_1_file)
    pair_2_df = import_fastq(pair_2_file)
    
    pair_2_df["reverse"] = pair_2_df["sequences"].apply(lambda x: reverse_complement_seq(x))
    pair_1_df["merged"] = pair_1_df["sequences"] + "NNNNNNNNNN" + pair_2_df["reverse"]
    pair_1_df["qual_merge"] = pair_1_df["quality"] + "##########" + pair_2_df["quality"].apply(lambda x: reverse_quality(x))
    
    final_df = pair_1_df[["ID", "merged", "junk", "qual_merge"]]
    final_df["ID"] = "@" + final_df["ID"]
    
    final_df.to_csv(export_file, index = False, header = False, quoting = 3, sep = "\n")
    
    