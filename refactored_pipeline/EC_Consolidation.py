#!/usr/bin/env python

import os
import os.path
import sys
import shutil

Input_File = sys.argv[1]
Input_Name = os.path.splitext(os.path.basename(Input_File))[0]
Output_Dir = sys.argv[2]

SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"
SWISS_PROT_MAP = "/home/j/jparkins/mobolaji/Databases/SwissProt_EC_Mapping.tsv"

detect_dir = os.path.join(Output_Dir, "Detect")
priam_dir = os.path.join(Output_Dir, "PRIAM")
diamond_dir = os.path.join(Output_Dir, "Diamond")
con_dir = os.path.join(Output_Dir, "Consolidated")
try:
    os.mkdir(con_dir)
except:
    shutil.rmtree(con_dir)
    os.mkdir(con_dir)

mapping_dict = {}
with open(SWISS_PROT_MAP, "r") as mapping:
    for line in mapping.readlines():
        line_as_list = line.split("\t")
        mapping_dict[line_as_list[0]] = set(line_as_list[2:])

detect_ECs = os.path.join(detect_dir, Input_Name + ".toppred")
with open(detect_ECs, "r") as topred:
    with open(os.path.join(detect_dir, Input_Name + ".toppred.cutoff"), "w") as cutoff:
        for line in topred.readlines():
            line_as_list = line.split("\t")
            if line_as_list[2] == "probability":
                continue
            if float(line_as_list[2]) >= 0.2 and int(line_as_list[3]) > 5:
                cutoff.write(line)

priam_ECs = os.path.join(priam_dir, Input_Name + ".ECs")
with open(os.path.join(priam_dir, "RESULTS", "paj_" + Input_Name.split("_proteins")[0] + "_PRIAM" + "_seqsECs.tab"), "r") as ECs:
    with open(priam_ECs, "w") as processedECs:
        for line in ECs.readlines():
            line_as_list = line.split("\t")
            if len(line_as_list[1].split(";")) > 1:
                for EC in range(len(line_as_list[1].split(";"))):
                    processedECs.write("\t".join([line_as_list[0], line_as_list[1].split(";")[EC].strip()]))
                    processedECs.write("\n")
            else:
                processedECs.write(line)

diamond_ECs = os.path.join(diamond_dir, Input_Name + ".ECs")
with open(os.path.join(diamond_dir, Input_Name.split("_proteins")[0] + ".blastout"), "r") as blastout:
    with open(diamond_ECs, "w") as ecout:
        for line in blastout.readlines():
            line_as_list = line.strip().split("\t")
            for EC in mapping_dict:
                if line_as_list[1] in mapping_dict[EC]:
                    ecout.write("\t".join([line_as_list[0], EC + "\n"]))

with open(os.path.join(con_dir, Input_Name + ".ECs_PB"), "w") as PB_out:
    with open(priam_ECs, "r") as priam_ECs_in:
        priam_preds = priam_ECs_in.readlines()
    with open(diamond_ECs, "r") as diamond_ECs_in:
        diamond_preds = diamond_ECs_in.readlines()
    PB_preds = []
    for priam_ec in priam_preds:
        if priam_ec in diamond_preds:
            PB_preds.append(priam_ec)
    PB_out.writelines(PB_preds)
with open(detect_ECs, "r") as detect_ECs_in:
    detect_preds = []
    for line in detect_ECs_in.readlines():
        line_as_list = line.split("\t")
        line_as_list = "\t".join(line_as_list[:2]) + "\n"
        detect_preds.append(line_as_list)
All_preds = set()
for pred in detect_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
for pred in PB_preds:
    if len(pred.split("\t")[1].split(".")) == 4:
        All_preds.add(pred)
with open(os.path.join(con_dir, Input_Name + ".ECs_All"), "w") as ec_out:
    ec_out.writelines(sorted(All_preds))