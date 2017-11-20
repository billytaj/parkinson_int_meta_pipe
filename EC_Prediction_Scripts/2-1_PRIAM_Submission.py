#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess

Input_Dir = sys.argv[1]
Output_Dir = sys.argv[2]

PRIAM = "/home/j/jparkins/mobolaji/Tools/PRIAM/PRIAM_search.jar"
BLAST = "/home/j/jparkins/mobolaji/Tools/BLAST/blast-2.2.26/bin/"
Python = "/home/j/jparkins/mobolaji/python"

PRIAM_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N NAME
cd OUT_PRIAM
module load java
java -jar PRIAM_PROG -n NAME -i INPUT -p PRIAM_PROF -od OUT_PRIAM -e T -pt 0.5 -mo -1 -mp 70 -cc T -cg T -bd BLAST"""

Watcher_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N PRIAM_Watcher
cd $PBS_O_WORKDIR
PYTHON WATCHER INPUT OUTPUT PYTHON"""

priam_dir_top = os.path.join(Output_Dir, "priam")
try:
    os.mkdir(priam_dir_top)
except:
    shutil.rmtree(priam_dir_top)
    os.mkdir(priam_dir_top)

for infile in os.listdir(Input_Dir):
    fname = os.path.splitext(infile)[0]
    priam_dir = os.path.join(priam_dir_top, fname)
    try:
        os.mkdir(priam_dir)
    except:
        shutil.rmtree(priam_dir)
        os.mkdir(priam_dir)
    with open(os.path.join(priam_dir, fname + ".pbs"), "w") as priam_script_out:
        for line in PRIAM_Script.splitlines():
            if "NAME" in line:
                line = line.replace("NAME", fname + "_PRIAM")
            if "PRIAM_PROG" in line:
                line = line.replace("PRIAM_PROG", PRIAM)
            if "BLAST" in line:
                line = line.replace("BLAST", BLAST)
            if "INPUT" in line:
                line = line.replace("INPUT", os.path.join(Input_Dir, infile))
            if "PRIAM_PROF" in line:
                line = line.replace("PRIAM_PROF", os.path.join(os.path.dirname(PRIAM), "PRIAM_MAR15"))
            if "OUT_PRIAM" in line:
                line = line.replace("OUT_PRIAM", priam_dir)
            priam_script_out.write(line + "\n")
    subprocess.call(["qsub", os.path.join(priam_dir, os.path.splitext(infile)[0] + ".pbs")])


with open(os.path.join(Output_Dir, "PRIAM_Watcher.pbs"), "w") as watcher_script_out:
    Watcher = os.path.join(os.path.dirname(sys.argv[0]), "2-2_PRIAM_Postprocessing.py")
    for line in Watcher_Script.splitlines():
        if "PYTHON" in line:
            line = line.replace("PYTHON", Python)
        if "WATCHER" in line:
            line = line.replace("WATCHER", Watcher)
        if "INPUT" in line:
            line = line.replace("INPUT", Input_Dir)
        if "OUTPUT" in line:
            line = line.replace("OUTPUT", Output_Dir)
        watcher_script_out.write(line + "\n")
subprocess.call(["qsub", os.path.join(Output_Dir, "PRIAM_Watcher.pbs")])