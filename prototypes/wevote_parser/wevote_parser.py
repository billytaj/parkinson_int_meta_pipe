import os
import sys
import pandas as pd

if __name__ == "__main__":
    wevote_df = pd.read_csv(sys.argv[1], sep = "\t", header = None)
    #wevote_df.dropna(inplace = True)
    wevote_df.columns = ["read_id", "tool_count", "tools_can_classify", "agreed", "score", "junk", "tool_0", "tool_1", "tool_2", "tool_3", "tool_4", "assignment"]
    
    print(wevote_df)
    wevote_df = wevote_df[["read_id", "assignment"]]
    print(wevote_df)
    
    wevote_df.to_csv("wevote_results.tsv", sep = "\t", index = False, header = None)