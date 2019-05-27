#this collects the @SQ portion of the samfile that we care about.  we want read lengths

import os
import sys
import pandas as pd

if __name__ == "__main__":

    
    sam_file = open(sys.argv[1], "r")
    line_count = 0
    list_col = []
    for line in sam_file:
        if (line.startswith("@SQ")):
            
            
            line_ID = line.split("\t")[1].split("SN:")[1].split("|")
            #print(line_ID)
            read_length = line_ID[-1]
            line_ID.pop(-1)
            new_ID = ""
            for item in line_ID:
                new_ID += item + "|"
            new_ID = new_ID[:-1]
            #print("new ID:", new_ID)
            #print("read length:", read_length)
            line_count += 1
            data_list = [new_ID, read_length]
            list_col.append(data_list)
        elif(line.startswith("@PQ") or line.startswith("@RQ")):
            break
    #print("line count:", line_count)
    
    read_length_df = pd.DataFrame(list_col)
    read_length_df.columns = ["GeneID", "Length"]
    new_export_name = os.path.splitext(sys.argv[1])[0] + "_read_length.csv"
    print(new_export_name)
    
    read_length_df.to_csv(new_export_name, mode = "w", index = False)
    #print(read_length_df)
    
    
    
