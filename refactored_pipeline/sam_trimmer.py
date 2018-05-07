import os
import sys
# this module strips away the pesky optional segments of the SAM file format.
# this is used so that pandas' read_csv importing is clean

if __name__ == "__main__":
    sam_path = sys.argv[1]
    out_path = sys.argv[2]
    
    sam_file = open(sam_path, "r")
    out_file = open(out_path, "w")
    count = 0
    for line in sam_file:
        cur_line = line.split('\t')
        line_length = len(cur_line)
        if(line_length > 12):
            print("line[", count, "]: length:", line_length)
            print(cur_line)
            print("--------------------------------------------------")
        count += 1
        if(line.startswith("@") or len(line) <= 1 or line_length < 13):
            continue
        new_line = cur_line[0:12]
        new_line[11] = new_line[11]+'\n'
        new_line = '\t'.join(new_line)
        out_file.write(new_line)