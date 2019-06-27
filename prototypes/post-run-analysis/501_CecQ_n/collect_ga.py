import pandas as pd
import os
import sys

if __name__ == "__main__":
    report_df = pd.read_csv(sys.argv[1], sep = "\t", error_bad_lines = False)
    
    ID_df = report_df["GeneID"].str.split("|", expand = True)
    ID_df["reads"] = report_df["Reads"]
    print(ID_df)
    #6 = the partial taxa
    group_df = ID_df[6].to_frame()
    group_df["reads"] = report_df["Reads"]
    group_df.columns = (["taxa", "reads"])
    print(group_df)
    group_df = group_df.groupby(["taxa"], as_index = False).sum()
    print("----------------------------")
    print(group_df[group_df["taxa"] == "g__Parabacteroides.s__Parabacteroides_sp_ASF519"])
    