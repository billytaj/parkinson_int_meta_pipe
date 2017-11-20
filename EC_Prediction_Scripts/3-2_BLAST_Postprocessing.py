#!/usr/bin/env python

import os
import os.path
import sys
import shutil
import time

Input_Dir = sys.argv[1]
Output_Dir = sys.argv[2]
Python = sys.argv[3]
SWISS_PROT_MAP = sys.argv[4]

done = False

while done is False:
    for file in os.listdir(os.getcwd()):
        if file.startswith("BLAST-DIAMOND"):
            done = True
            break
    else:
        time.sleep(30)

for file in os.listdir(os.getcwd()):
    if file.startswith("BLAST-DIAMOND"):
        try:
            os.rename(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "blast", file))
        except:
            shutil.copy(os.path.join(os.getcwd(), file), os.path.join(Output_Dir, "blast", file))

mapping_dict = {}

with open(SWISS_PROT_MAP, "r") as mapping:
    for line in mapping.readlines():
        line_as_list = line.split("\t")
        mapping_dict[line_as_list[0]] = set(line_as_list[2:])

file_list = []

for infile in os.listdir(Input_Dir):
    fname = os.path.splitext(infile)[0]
    file_list.append(fname)

for fname in file_list:
    with open(os.path.join(os.path.join(Output_Dir, "blast", fname, fname + ".blastout")), "r") as blastout:
        with open(os.path.join(os.path.join(Output_Dir, "blast", fname, fname + ".ECs")), "w") as ecout:
            for line in blastout.readlines():
                line_as_list = line.strip().split("\t")
                for EC in mapping_dict:
                    if line_as_list[1] in mapping_dict[EC]:
                        ecout.write("\t".join([line_as_list[0], EC + "\n"]))