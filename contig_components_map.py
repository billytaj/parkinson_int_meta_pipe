
import sys
import pandas as pd
pair_1_sam_path = sys.argv[1]

pair_1_sam_df = pd.read_csv(pair_1_sam_path, error_bad_lines=False, header=None, sep="\t")
pair_1_sam_df.iloc[:, 1] = pair_1_sam_df.iloc[:, 1].apply(lambda x: bin(int(x))[2:].zfill(11)[8])

selected_mapped = pair_1_sam_df.loc[pair_1_sam_df.iloc[:, 1] == "0"].iloc[:, 0:3]
selected_mapped.columns = ["reads", "flag", "contig"]


selected_mapped.index = selected_mapped.groupby('contig').cumcount() # this looks like it makes the indices, by finding the longest row
#print(selected_mapped.pivot(values='reads', columns='contig'))

t2 = selected_mapped.pivot(values='reads', columns='contig')

new_t2 = t2.T
new_t2["freq"] = new_t2.apply(lambda x: x.count(), axis=1)
cols = new_t2.columns.tolist()
cols = cols[-1:] + cols[:-1]
new_t2 = new_t2[cols]
new_t2.to_csv("contig_map.tsv", sep = '\t', mode = "w+", header=False)