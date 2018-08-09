import pandas as pd
import sys
import os
import numpy as np

def string_to_ascii_array(line):
    new_line = ""
    for item in line:
        new_line += str(ord(item)) + ","
    return new_line[:-1]

if __name__ == "__main__":
    df_file = pd.read_csv(sys.argv[1], header = None, names=None, sep='\n', skip_blank_lines=False)
    df_0 = pd.DataFrame(df_file.values.reshape(int(len(df_file)/4), 4))
    df_0.columns = ["ID", "seq", "junk", "quality"]
    
    df_0.drop(["ID", "seq", "junk"], axis =1, inplace=True)
    df_0["quality"] = df_0["quality"].apply(lambda x: string_to_ascii_array(x))
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
    new_name = os.path.split(sys.argv[1])[1].split(".")[0] + "_base_pair_quality_report.csv"
    stats_df.to_csv(new_name, mode="w+")
    #df_0.boxplot(figsize = [80,20])
    #print(stats_df)
    #print(df_0)
    
    