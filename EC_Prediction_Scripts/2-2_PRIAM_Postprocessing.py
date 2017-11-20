#!/usr/bin/env python

import os
import os.path
import sys
import shutil
import time

Input_Dir = sys.argv[1]
Output_Dir = sys.argv[2]
Python = sys.argv[3]

done = False

file_list = []

for infile in os.listdir(Input_Dir):
    fname = os.path.splitext(infile)[0]
    file_list.append(fname)

while done is False:
    for fname in file_list:
        for file in os.listdir(os.getcwd()):
            if file.startswith(fname + "_PRIAM"):
                break
        else:
            time.sleep(30)
            break
    else:
        done = True

for fname in file_list:
    for file in os.listdir(os.getcwd()):
        if file.startswith(fname + "_PRIAM"):
            try:
                os.rename(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "priam", fname, file))
            except:
                shutil.copy(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "priam", fname, file))
    with open(os.path.join(Output_Dir, "priam", fname, "RESULTS", "paj_" + fname + "_PRIAM" + "_seqsECs.tab"), "r") as ECs:
        with open(os.path.join(Output_Dir, "priam", fname, fname + ".ECs"), "w") as processedECs:
            for line in ECs.readlines():
                line_as_list = line.split("\t")
                if len(line_as_list[1].split(";")) > 1:
                    for EC in range(len(line_as_list[1].split(";"))):
                        processedECs.write("\t".join([line_as_list[0], line_as_list[1].split(";")[EC].strip()]))
                        processedECs.write("\n")
                else:
                    processedECs.write(line)