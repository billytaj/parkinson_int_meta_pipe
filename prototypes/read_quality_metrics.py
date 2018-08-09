import pandas as pd
import sys
import os
import numpy as np

class read_quality_metrics:
    def __init__(self, fastq_file):
        self.filename = fastq_file
        self.df_file = pd.read_csv(fastq_file, header = None, names=None, sep='\n', skip_blank_lines=False)
        self.df_orig = pd.DataFrame(self.df_file.values.reshape(int(len(self.df_file)/4), 4))
        self.df_orig.columns = ["ID", "seq", "junk", "quality"]
        
    def string_to_ascii_array(self, line):
        new_line = ""
        for item in line:
            new_line += str(ord(item)) + ","
        return new_line[:-1]
    
    def avg_ascii(self, line):
        sum = 0
        for item in line:
            sum += ord(item)
        return sum / len(line)
        
    
    
    def per_base_quality(self):
        df_0 = self.df_orig["quality"]
        df_0["quality"] = df_0["quality"].apply(lambda x: self.string_to_ascii_array(x))
        df_0 = df_0["quality"].str.split(",", expand=True).rename(columns = lambda x: "bp_" + str(x))
        df_0 = df_0.apply(pd.to_numeric)
        df_0.loc["avg"] = df_0.select_dtypes(pd.np.number).sum() / df_0.shape[0] #doing on select_dtypes means we only consider numerics.  
        
        stats_df = pd.DataFrame(df_0.loc["avg"]).transpose()
        stats_df.loc["SD"] = df_0.select_dtypes(pd.np.number).std()
        stats_df.loc["Median"] = df_0.select_dtypes(pd.np.number).median()
        stats_df.loc["Max"] = df_0.select_dtypes(pd.np.number).max()
        stats_df.loc["Min"] = df_0.select_dtypes(pd.np.number).min()
        stats_df.loc["Q1"] = df_0.select_dtypes(pd.np.number).quantile(0.25)
        stats_df.loc["Q3"] = df_0.select_dtypes(pd.np.number).quantile(0.75)
        new_name = os.path.split(self.filename)[1].split(".")[0] + "_per_base_quality_report.csv"
        stats_df.to_csv(new_name, mode="w+", header=False, index=False)
        
    def per_sequence_quality(self):
        #number of reads vs mean quality of base pairs in the read
        df_0 = pd.DataFrame(self.df_orig["quality"])
        df_0["len"] = df_0["quality"].apply(lambda x: len(x))
        df_0["quality"] = df_0["quality"].apply(lambda x: self.avg_ascii(x))
        df_0.hist(column="quality")#, by = "len")
        new_name = os.path.split(self.filename)[1].split(".")[0] + "_per_seq_quality_report.csv"
        df_0.to_csv(new_name, mode = "w+", header=False, index=False)
        
        
if __name__ == "__main__":
    fastq_file = sys.argv[1]
    read_stats_obj = read_quality_metrics(fastq_file)
    
    read_stats_obj.per_sequence_quality()
    