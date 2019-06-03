import sys
import os
import pandas as pd


def import_names_file(file):
    names_df = pd.read_csv(file, sep = "\t", error_bad_lines = False)
    names_df = names_df.groupby(["English name"], as_index = False).sum()
    new_name = file.split(".tsv")[0] + "_summary.csv"
    names_df.to_csv(new_name, index = False)
    return names_df


if __name__ == "__main__":
    names_file = sys.argv[1]
    
    names_df = import_names_file(names_file)
    
    print(names_df)
    