import os
import pandas as pd
import sys
import subprocess as sp

def extract_rRNA_ID(inf_file):
    #If there's a better way to do this, I'd like to see it.
    ID_list = []
    inf_list = open(inf_file, mode='r')
    
    for item in inf_list:
        if(not item.startswith("#")):
            for subitem in item.split(' '):
                if(subitem.startswith('ERR')):
                    ID_list.append(subitem)
                    break
                
        elif(len(item) == 2):
            break
    return ID_list

def filter_rRNA(rRNA_ID_list, fastq_sequence):
    fastq_df = pd.read_csv(fastq_sequence, header=None, names = [None])
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    
    
if __name__ == "__main__":
    print("hey")
    
    inf_file = sys.argv[1]
    fastq_sequence = sys.argv[2]
    ID_list = extract_rRNA_ID(inf_file)
    
    filter_rRNA(ID_list, fastq_sequence)
    
    
    """
    p = sp.Popen(inf_file, stdout = sp.PIPE, shell=True)
    print(p)
    """
    
    """
    spec_read = pd.read_csv(p.stdout, usecols=[
    '1',
    '2',
    '3',
    '4',
    '5',
    '6',
    '7',
    '8',
    '9',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18'])
    """