import pandas as pd
import sys
import time
import numpy as np

def sort_and_export(input_file, output_file):

    start_total_call = time.clock()
    print("input file:", input_file)
    df = pd.read_csv(input_file, header=None, names=[None], sep='\n', skip_blank_lines = False)
    end_read_time = time.clock()
    
    # reshaping, aka: unflattening the 1-D array into something meaningful to us
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))

    df.columns = ["ID", "sequences", "junk", "quality"]
    #df["sort_ID"] = df["ID"]
    #df["sort_ID"] = df["sort_ID"].str.replace('@ERR', '')
    # can be done in a single line.  
    # what's happening: we're extracting the read number from the ID column, and making a new column out of it.
    #df["sort_ID"] = df["sort_ID"].str.split(' ', n=1, expand=True)[0]
    #df["sort_ID"] = df["sort_ID"].str.split('.', n=1, expand=True)[1].astype('int')
    # next, we'll sort the dataframe with the extract read number
    #df = df.sort_values(by=['sort_ID'])
    # since we don't need it in the final output, we're erasing the column
    #df = df.drop(['sort_ID'], axis=1)
    
    # old code just strips away the original ID, and writes that 
    #df["ID"] = df["ID"].str.split(' ', n=1, expand=True)[0]
    df = df.sort_values(by=['ID'])
    end_df_time = time.clock()
    df.to_csv(output_file, sep='\n', mode = 'w+', header=False, index=False)
    end_total_call = time.clock()
    #-------------------------------------------------------------------------------------------
    
    print("import csv time:", end_read_time - start_total_call, "s")
    print("dataframe interpret time:", end_df_time - end_read_time, "s")
    print("write time:", end_total_call - end_df_time, "s")
    print("total runtime:", end_total_call - start_total_call, "s")

    
if __name__ == "__main__":
    scinet_flag = False
    if(len(sys.argv) < 3):
        print("something wrong with the args")
        sys.exit()
    else:
        input_fastq_path = sys.argv[1]
        export_path = sys.argv[2]
    
    sort_and_export(input_fastq_path, export_path)
