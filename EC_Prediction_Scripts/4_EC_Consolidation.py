#!/usr/bin/env python

import os
import os.path
import sys
import shutil

#ADD WATCHER

Input_Dir = sys.argv[1]
Split_Dir = sys.argv[2]
Output_Dir = sys.argv[3]

SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"
SWISS_PROT_MAP = "/home/j/jparkins/mobolaji/Databases/SwissProt_EC_Mapping.tsv"
BLAST_Plus = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/"
DIAMOND = "/home/j/jparkins/mobolaji/Tools/Diamond/diamond"
Python = "/home/j/jparkins/mobolaji/python"

detect_dir_top = os.path.join(Output_Dir, "detect")
priam_dir_top = os.path.join(Output_Dir, "priam")
blast_dir_top = os.path.join(Output_Dir, "blast")
con_dir_top = os.path.join(Output_Dir, "consolidated")
try:
    os.mkdir(con_dir_top)
except:
    shutil.rmtree(con_dir_top)
    os.mkdir(con_dir_top)

for infile in os.listdir(Split_Dir):
    fname = os.path.splitext(infile)[0]
    detect_dir = os.path.join(detect_dir_top, fname)
    detect_ECs = os.path.join(detect_dir, fname + ".topred.cutoff")
    priam_dir = os.path.join(priam_dir_top, fname)
    priam_ECs = os.path.join(priam_dir, fname + ".ECs")
    blast_dir = os.path.join(blast_dir_top, fname)
    blast_ECs = os.path.join(blast_dir, fname + ".ECs")
    con_dir = os.path.join(con_dir_top, fname)
    try:
        os.mkdir(con_dir)
    except:
        shutil.rmtree(con_dir)
        os.mkdir(con_dir)
    with open(os.path.join(con_dir, fname + ".ECs_PB"), "w") as PB_out:
        with open(priam_ECs, "r") as priam_ECs_in:
            priam_preds = priam_ECs_in.readlines()
        with open(blast_ECs, "r") as blast_ECs_in:
            blast_preds = blast_ECs_in.readlines()
        PB_preds = []
        for priam_ec in priam_preds:
            if priam_ec in blast_preds:
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
    with open(os.path.join(con_dir, fname + ".ECs_All"), "w") as ec_out:
        ec_out.writelines(sorted(All_preds))

for infile in os.listdir(Input_Dir):
	if os.path.isdir(os.path.join(Input_Dir, infile)):
		continue
    fname = os.path.splitext(infile)[0]
    Out = []
    for folder in os.listdir(con_dir_top):
        if folder.startswith(fname):
            with open(os.path.join(con_dir_top, folder, folder + ".ECs_All"), "r") as EC_Split:
                Out.extend(EC_Split.readlines())
    else:
        for index, line in enumerate(Out):
            if not line.endswith("\n"):
                Out[index] = line + "\n"
        with open(os.path.join(con_dir_top, fname + ".ECs_All"), "w") as ECs:
                ECs.writelines(sorted(Out))