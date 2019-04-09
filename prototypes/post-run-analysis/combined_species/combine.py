import os
import sys
import pandas as pd

#this code needs to merge all the results from the reads together, to get a final table.
#We'll do it be the average

#It's gonna be a concat + groupby.mean()

if __name__ == "__main__":
    input_dir = sys.argv[1]
    output_dir = sys.argv[1]
    list_of_files = os.listdir(input_dir)
    
    count = 0
    master_df = None
    for item in list_of_files:
        if(item.endswith(".csv")):
            print("working on:", item)
            if(count == 0):
                count += 1
                master_df = pd.read_csv(item, error_bad_lines = False)
                
            else:
                new_df = pd.read_csv(item, error_bad_lines = False)
                master_df = pd.concat((master_df, new_df))
                
    master_df = master_df.groupby(["taxa"], as_index = False).mean()
    print(master_df)
    new_name = os.path.join(output_dir, "combined_avg.tsv")
    master_df.to_csv(new_name, sep = "\t", index = False)
    