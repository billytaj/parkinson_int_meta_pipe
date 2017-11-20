#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess

Input_Dir = sys.argv[1]
Output_Dir = sys.argv[2]

DETECT = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
Python = "/home/j/jparkins/mobolaji/python"

Detect_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N NAME
cd OUT0
export PATH=/home/j/jparkins/ctorma/emboss/bin/:$PATH
module load gcc blast extras
module load intel/15.0.2
module load python/2.7.8
PYTHON DETECT INPUT --output_file OUT1 --top_predictions_file OUT2 --num_threads 8"""

Watcher_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N Detect_Watcher
cd $PBS_O_WORKDIR
PYTHON WATCHER INPUT OUTPUT PYTHON"""

detect_dir_top = os.path.join(Output_Dir, "detect")
try:
    os.mkdir(detect_dir_top)
except:
    shutil.rmtree(detect_dir_top)
    os.mkdir(detect_dir_top)

for infile in os.listdir(Input_Dir):
    fname = os.path.splitext(infile)[0]
    detect_dir = os.path.join(detect_dir_top, fname)
    try:
        os.mkdir(detect_dir)
    except:
        shutil.rmtree(detect_dir)
        os.mkdir(detect_dir)
    with open(os.path.join(detect_dir, fname + ".pbs"), "w") as detect_script_out:
        for line in Detect_Script.splitlines():
            if "NAME" in line:
                line = line.replace("NAME", fname +"_Detect")
            if "PYTHON" in line:
                line = line.replace("PYTHON", Python)
            if "DETECT" in line:
                line = line.replace("DETECT", DETECT)
            if "INPUT" in line:
                line = line.replace("INPUT", os.path.join(Input_Dir, infile))
            if "OUT0" in line:
                line = line.replace("OUT0", detect_dir)
            if "OUT1" in line:
                line = line.replace("OUT1", os.path.join(detect_dir, fname + ".detect"))
            if "OUT2" in line:
                line = line.replace("OUT2", os.path.join(detect_dir, fname + ".topred"))
            detect_script_out.write(line + "\n")
    subprocess.call(["qsub", os.path.join(detect_dir, os.path.splitext(infile)[0] + ".pbs")])

with open(os.path.join(Output_Dir, "Detect_Watcher.pbs"), "w") as watcher_script_out:
    Watcher = os.path.join(os.path.dirname(sys.argv[0]), "1-2_Detect_Postprocessing.py")
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
subprocess.call(["qsub", os.path.join(Output_Dir, "Detect_Watcher.pbs")])