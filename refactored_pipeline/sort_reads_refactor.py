import pandas as pd
import sys
import time
import numpy as np


input_fastq = sys.argv[1]

start_total_call = time.clock()
df = pd.read_csv(input_fastq, header=None, names=[None])
end_read_time = time.clock()
print("import csv time:", end_read_time - start_total_call, "s")
#print("length:", len(df))
df = df.values.reshape(int(len(df)/4), 4)
#print(full_df[0:10])
new_df = pd.DataFrame(df)
new_df.columns = ["ID", "sequences", "junk", "quality"]
#new_df = new_df.drop(['junk'], axis = 1)
new_df = new_df.sort_values(by=['ID'])
end_df_time = time.clock()
print("dataframe interpret time:", end_df_time - end_read_time, "s")
#print(new_df)
export_filename = "sorted_" + input_fastq
new_df.to_csv(export_filename, sep='\n', mode = 'w+', header=False, index=False)
end_total_call = time.clock()
print("total runtime:", end_total_call - start_total_call, "s")

