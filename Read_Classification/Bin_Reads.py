#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
import multiprocessing
from Bio import SeqIO

Python = "/home/j/jparkins/mobolaji/python"
BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
Sort_Reads = "/home/j/jparkins/mobolaji/Read_Classification/Sort_Reads.py"

Kaiju = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju"
Centrifuge = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge"
Classification_postprocess = "/home/j/jparkins/mobolaji/Read_Classification/Lineage_consensus.py"

Metatranscriptome_Pipeline = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Metatranscriptome_pipeline.py"

Threads = str(multiprocessing.cpu_count())

PBS_Submit = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N NAME

module load gcc intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

input_folder = sys.argv[1]
output_folder = sys.argv[2]
classification_rank = sys.argv[3]
cutoff = sys.argv[4]
try:
    Test_run = sys.argv[5]
except:
    Test_run = ""


for genome in sorted(os.listdir(input_folder)):
    if genome.endswith("1.fastq"):
        genome_name = genome[:-7]
        Input_File = os.path.join(output_folder, genome_name + ".fastq")
        Input_Filepath = os.path.splitext(Input_File)[0]
        Input_File1 = Input_Filepath + "1"
        Input_File2 = Input_Filepath + "2"
        Input_Path = os.path.dirname(Input_File)
        Input_FName = os.path.basename(Input_File)
        try:
            File_stats = subprocess.check_output([os.path.join(BBMap_Dir, "testformat.sh"), os.path.join(input_folder, genome_name + "1.fastq")])
        except:
            sys.exit("Use module load java before running script")
        if File_stats.startswith("sanger"):
            Qual = "33"
        else:
            Qual = "64"

        os.chdir(Input_Path)

        # Taxonomic Annotation
        COMMANDS = [
        Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db.fmi" + " -i " + Input_File1 + ".fastq" + " -j " + Input_File2 + ".fastq" + " -z " + Threads + " -o " + Input_Filepath + "_KaijuOut.tsv",
        Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/p+h+v" + " -1 " + Input_File1 + ".fastq" + " -2 " + Input_File2 + ".fastq" + " --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_CentrifugeOut.tsv" + " --report-file " + "/dev/null",
        "mkdir -p " + Input_Filepath + "_Bins",
        Python + " " + Classification_postprocess + " " + classification_rank + " " + cutoff + " " + Input_File1 + ".fastq" + " " + Input_File2 + ".fastq" + " " + Input_Filepath + "_Bins" + " " + Input_Filepath + "_read_classification.tsv" + " " + Input_Filepath + "_read_classification_summary.tsv" + " " + Input_Filepath + "_KaijuOut.tsv" + " " + Input_Filepath + "_CentrifugeOut.tsv"
        ]

        if Test_run != "":
            COMMANDS.extend(["ssh gpc01" + " " + "\"" + "module load java;" + " " + Python + " " + Metatranscriptome_Pipeline + " " + Input_Filepath + "_Bins" + " " + output_folder + "\""])

        with open(os.path.splitext(Input_FName)[0] + "_Classify.pbs", "w") as PBS_script_out:
            for line in PBS_Submit.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Classify")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDS))
                    break
                PBS_script_out.write(line + "\n")
        JobID = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Classify.pbs"])
