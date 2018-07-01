import sys
import pandas as pd
fasta_path = sys.argv[1]

fasta_df = pd.read_csv(fasta_path, error_bad_lines=False, header=None, sep="\n")                #import the fasta
fasta_df.columns = ["row"]                                                                              
new_df = pd.DataFrame(fasta_df.loc[fasta_df.row.str.contains('>')])                             #grab all the IDs
new_df.columns =["names"]
new_data_df = fasta_df.loc[~fasta_df.row.str.contains('>')]                                     #grab the data
new_data_df.columns = ["data"]
fasta_df = new_df.join(new_data_df, how='outer')                                                #join them into a 2-col DF
print(fasta_df.info())
fasta_df["names"] = fasta_df.fillna(method='ffill')                                             #fill in the blank spaces in the name section
fasta_df.dropna(inplace=True)                                                                   #remove all rows with no sequences
fasta_df.index = fasta_df.groupby('names').cumcount()                                           #index it for transform
temp_columns = fasta_df.index                                                                   #save the index names for later
fasta_df = fasta_df.pivot(values='data', columns = 'names')                                     #pivot
fasta_df = fasta_df.T                                                                           #transpose
fasta_df["sequence"] = fasta_df[fasta_df.columns[:]].apply(lambda x: "".join(x.dropna()),axis=1)#consolidate all cols into a single sequence
fasta_df.drop(temp_columns, axis=1, inplace=True)                                               #drop all other columns
fasta_df.to_csv("out_new.csv", sep=",", mode="w+")                                              #export (testing)