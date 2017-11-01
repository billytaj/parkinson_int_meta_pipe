import pandas as pd
import sys
import time
import numpy as np

def main(input_file, scinet_flag):

    export_filepath = ""
    if(scinet_flag):
        export_filepath = "/scratch/j/jparkins/billyc59/sorted_" +  input_file.split('/')[-1]
        
    elif "\\" in input_file:
        split_filepath = input_file.split('\\') #assumes 1 layer down
        export_filepath = split_filepath[0] + "\sorted_" + split_filepath[1]
    else:
        export_filepath = "sorted_" + input_file
    
    start_total_call = time.clock()
    df = pd.read_csv(input_file, header=None, names=[None])
    end_read_time = time.clock()
    
    # reshaping, aka: unflattening the 1-D array into something meaningful to us
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))

    df.columns = ["ID", "sequences", "junk", "quality"]
    df["sort_ID"] = df["ID"]
    df["sort_ID"] = df["sort_ID"].str.replace('@ERR', '')
    # can be done in a single line.  
    # what's happening: we're extracting the read number from the ID column, and making a new column out of it.
    df["sort_ID"] = df["sort_ID"].str.split(' ', n=1, expand=True)[0]
    df["sort_ID"] = df["sort_ID"].str.split('.', n=1, expand=True)[1].astype('int')
    # next, we'll sort the dataframe with the extract read number
    df = df.sort_values(by=['sort_ID'])
    # since we don't need it in the final output, we're erasing the column
    df = df.drop(['sort_ID'], axis=1)
    end_df_time = time.clock()
    df.to_csv(export_filepath, sep='\n', mode = 'w+', header=False, index=False)
    end_total_call = time.clock()
    #-------------------------------------------------------------------------------------------
    
    print("import csv time:", end_read_time - start_total_call, "s")
    print("dataframe interpret time:", end_df_time - end_read_time, "s")
    print("write time:", end_total_call - end_df_time, "s")
    print("total runtime:", end_total_call - start_total_call, "s")

    
if __name__ == "__main__":
    scinet_flag = False
    if(len(sys.argv) < 2):
        print("missing launch mode flag (\"scinet\") or something else")
        sys.exit()
    else:
        input_fastq_path = sys.argv[1]
        launch_mode = sys.argv[2]
    if(launch_mode == "scinet"):
        scinet_flag = True
    main(input_fastq_path, scinet_flag)
