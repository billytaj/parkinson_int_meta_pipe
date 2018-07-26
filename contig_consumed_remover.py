#June 21, 2018:  This is a prototype.  it works logically, but it's not fully integrated to work with the pipeline

import sys
import pandas as pd


def make_map(sam_file):
    pair_1_sam_df = pd.read_csv(sam_file, error_bad_lines=False, header=None, sep="\t")                 # read it
    pair_1_sam_df.iloc[:, 1] = pair_1_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])   # rewrite the CIGAR flag
    selected_mapped = pair_1_sam_df.loc[pair_1_sam_df.iloc[:, 1] == "0"].iloc[:, 0:3]                   # select what we need from the SAM.  mapped reads only
    selected_mapped.columns = ["reads", "flag", "contig"]                                               # name the columns for easy manipulation 
    selected_mapped.index = selected_mapped.groupby('contig').cumcount()                                # this looks like it makes the indices, by finding the longest row
    t2 = selected_mapped.pivot(values='reads', columns='contig')                                        # pivot it so the contig is the columns
    new_t2 = t2.T                                                                                       # then transpose
    return new_t2

def final_concat(df_list):
    final_df = pd.concat(df_list, axis=1)                           #do initial concat
    new_cols = []                                                   #rename the entire df columns.  Not doing this will mess up the final result
    new_cols.extend(range(0, len(final_df.columns)))
    final_df.columns = new_cols
    final_df["freq"] = final_df.apply(lambda x: x.count(), axis=1)  #add the freq
    cols = final_df.columns.tolist()                
    cols = cols[-1:] + cols[:-1]                                    #make it the first column
    final_df = final_df[cols]
    return final_df                                                 
    
if __name__ == "__main__":
    sam_0 = sys.argv[1]
    df_list = []
    for i in range(0, 10):
        df_list.append(make_map(sam_0))
        
        
    final_df = final_concat(df_list)    
    final_df.to_csv("this_new.tsv", sep = '\t', mode = "w+")
   