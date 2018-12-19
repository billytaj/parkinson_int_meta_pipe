import sys
import time
import pandas as pd

def remove_tag(input_file, output_file):
    df = pd.read_csv(input_file, header=None, names=[None], sep='\n', skip_blank_lines = False, quoting=3)
    end_read_time = time.clock()
    
    df = pd.DataFrame(df.values.reshape(int(len(df)/4), 4))
    df.columns = ["ID", "sequences", "junk", "quality"]
    df["ID"] = df["ID"].apply(lambda x: x.split("/")[0])
    df.to_csv(output_file, sep='\n', mode = 'w+', header=False, index=False)
    
if __name__ == "__main__":
    scinet_flag = False
    if(len(sys.argv) != 3):
        print("something wrong with the args")
        sys.exit()
    else:
        input_fastq_path = sys.argv[1]      # expects a path to the fastq file, including the name
        export_path = sys.argv[2]           # expects a place to 
    
    remove_tag(input_fastq_path, export_path)
