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
            if file.startswith(fname + "_Detect"):
                break
        else:
            time.sleep(30)
            break
    else:
        done = True

for fname in file_list:
    for file in os.listdir(os.getcwd()):
        if file.startswith(fname + "_Detect"):
            try:
                os.rename(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "detect", fname, file))
            except:
                shutil.copy(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "detect", fname, file))
    with open(os.path.join(Output_Dir, "detect", fname, fname + ".topred"), "r") as topred:
        with open(os.path.join(Output_Dir, "detect", fname, fname + ".topred.cutoff"), "w") as cutoff:
            for line in topred.readlines():
                line_as_list = line.split("\t")
                if line_as_list[2] == "probability":
                    continue
                if float(line_as_list[2]) >= 0.9 and int(line_as_list[3]) > 5:
                    cutoff.write(line)