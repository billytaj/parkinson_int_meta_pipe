#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess

Input_Dir = sys.argv[1]
Output_Dir = sys.argv[2]

SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"
SWISS_PROT_MAP = "/home/j/jparkins/mobolaji/Databases/SwissProt_EC_Mapping.tsv"
BLAST_Plus = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/"
DIAMOND = "/home/j/jparkins/mobolaji/Tools/Diamond/diamond"
Python = "/home/j/jparkins/mobolaji/python"

MODE = "DIAMOND" # BLAST or DIAMOND

if MODE == "BLAST":
    BLASTP = os.path.join(BLAST_Plus, "blastp")
elif MODE == "DIAMOND":
    BLASTP = DIAMOND + " blastp"

BLAST_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N NAME
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
module load gnu-parallel/20140622
parallel -j 8 <<EOF
BLAST_CMD
EOF"""
BLAST_CMD = """  cd BLAST_DIR; BLASTP -query INPUT -db SWISS_PROT -outfmt "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident" -out OUTPUT -evalue 0.0000000001 -max_target_seqs 1"""

Watcher_Script = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N BLAST_Watcher
cd $PBS_O_WORKDIR
PYTHON WATCHER INPUT OUTPUT PYTHON SWISS_PROT_MAP"""

blast_dir_top = os.path.join(Output_Dir, "blast")
try:
    os.mkdir(blast_dir_top)
except:
    shutil.rmtree(blast_dir_top)
    os.mkdir(blast_dir_top)

with open(os.path.join(Output_Dir, "BLAST.pbs"), "w") as blast_script_out:
    for line in BLAST_Script.splitlines():
        if "NAME" in line:
            line = line.replace("NAME", "BLAST-DIAMOND")
        if "BLAST_CMD" in line:
            for infile in os.listdir(Input_Dir):
                fname = os.path.splitext(infile)[0]
                blast_dir = os.path.join(blast_dir_top, fname)
                input_file = os.path.join(Input_Dir, infile)
                out_file = os.path.join(blast_dir, fname + ".blastout")
                try:
                    os.mkdir(blast_dir)
                except:
                    shutil.rmtree(blast_dir)
                    os.mkdir(blast_dir)
                command = BLAST_CMD.replace("BLAST_DIR", blast_dir).replace("BLASTP", BLASTP).replace("INPUT", input_file).replace("SWISS_PROT", SWISS_PROT).replace("OUTPUT", out_file)
                if MODE == "DIAMOND":
                    command = command.replace(" -", " --").replace("\"", "").replace("max_target_seqs", "max-target-seqs")
                blast_script_out.write(command + "\n")
            continue
        blast_script_out.write(line + "\n")
subprocess.call(["qsub", os.path.join(Output_Dir, "BLAST.pbs")])

with open(os.path.join(Output_Dir, "BLAST_Watcher.pbs"), "w") as watcher_script_out:
    Watcher = os.path.join(os.path.dirname(sys.argv[0]), "3-2_BLAST_Postprocessing.py")
    for line in Watcher_Script.splitlines():
        if "PYTHON" in line:
            line = line.replace("PYTHON", Python)
        if "WATCHER" in line:
            line = line.replace("WATCHER", Watcher)
        if "INPUT" in line:
            line = line.replace("INPUT", Input_Dir)
        if "OUTPUT" in line:
            line = line.replace("OUTPUT", Output_Dir)
        if "SWISS_PROT_MAP" in line:
            line = line.replace("SWISS_PROT_MAP", SWISS_PROT_MAP)
        watcher_script_out.write(line + "\n")
subprocess.call(["qsub", os.path.join(Output_Dir, "BLAST_Watcher.pbs")])