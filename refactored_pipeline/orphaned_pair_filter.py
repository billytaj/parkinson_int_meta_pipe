import pandas as pd
import numpy as np
import sys
import os
# what this code does:
# takes in 2 pair fastq files, splits them up by 
# to scale this, break up the inputs.
# we expect this code to be called multiple times
# the naming convention will be based off of the inputs
# <pair_0>_unique.fastq
# <pair_1>_unique.fastq
# <pair_0>_<pair_1>_similar.fastq
os.chdir("C:\playground\\")
print("hely")
sample_0 = [1, 1, 1, 1, 3, 3, 3, 3, 2]
sample_1 = [1, 1, 1, 1, 2, 2, 2, 2, 3]

pre_df_0 = pd.read_csv("C:/playground/pair_1_clean.fastq", header=None, names=[None])
pre_df_1 = pd.read_csv("C:/playground/pair_2_clean.fastq", header=None, names=[None])
df_0 = pd.DataFrame(pre_df_0.values.reshape(int(len(pre_df_0)/4), 4))
df_1 = pd.DataFrame(pre_df_1.values.reshape(int(len(pre_df_1)/4), 4))
df_0.columns = ["ID", "seq", "junk", "quality"]
df_1.columns = ["ID", "seq", "junk", "quality"]
common = df_0.merge(df_1, on=["ID"])

df_0[df_0.ID.isin(common.ID)].to_csv("df_0_similar.fastq", sep = '\n', mode = 'w+', header = False, index = False)
df_1[df_1.ID.isin(common.ID)].to_csv("df_1_similar.fastq", sep = '\n', mode = 'w+', header = False, index = False)
df_0[~df_0.ID.isin(common.ID)].to_csv("unique.fastq", sep='\n', mode = 'w+', header=False, index = False)
df_1[~df_1.ID.isin(common.ID)].to_csv("unique.fastq", sep='\n', mode = 'a', header=False, index = False)
"""
if __name__ == "__main__":
    if(len(sys.argv) < 2):
        print("Too few input arguements.  Not filtering for orphans")
    else:
        pair_0_path = sys.argv[1]
        pair_1_path = sys.argv[2]
"""