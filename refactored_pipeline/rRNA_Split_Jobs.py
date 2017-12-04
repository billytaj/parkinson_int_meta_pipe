# This is the high-level mode that calls the infernal program

#!/usr/bin/env python

import sys
import os
import csv
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import clock as clock
from datetime import datetime as dt
import mt_pipe_paths as mpp

start_all = clock()
def 
if __name__ == "__main__":
    


#Input_File = sys.argv[1]
#Input_Path = os.path.dirname(Input_File)
#Input_FName = os.path.basename(Input_File)
#JobID1 = sys.argv[2]
#Filter_rRNA = sys.argv[3]

#Python = "/home/j/jparkins/mobolaji/python"
#Filter_rRNA = "/home/j/jparkins/billyc59/parkinson_int_meta_piperRNA_Filter.py"

#originally set to 6:00:00
PBS_Submit_LowMem = """#!/bin/bash
#PBS -l nodes=1:ppn=16,walltime=00:15:00 -q sandy
#PBS -N NAME
#PBS -e ERROR
#PBS -o OUTPUT

module load gcc intel/15.0.2 openmpi java blast extras anaconda3/4.0.0
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

Preprocess_jobs = []
export_filepath = ""
start_unpaired_infernal = clock()
#calls for all parts of a split orphans final

"""
for split in sorted(os.listdir(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants"))):
    
    if split.endswith(".fastq"):
        Split_File = os.path.splitext(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", split))[0]

        COMMANDSx = [mpp.Python + " " + Filter_rRNA + " " + Split_File + ".fastq" + " " + Split_File + "_mRNA.fastq" + " " + Split_File + "_rRNA.fastq"]

        with open(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs"), "w") as PBS_script_out:
            for line in PBS_Submit_LowMem.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(split)[0] + "_rRNA_Filter")
                if "ERROR" in line:
                    line = line.replace("ERROR", os.path.splitext(split)[0] + "_rRNA_Filter_ERR")
                if "OUTPUT" in line:
                    line = line.replace("OUTPUT", os.path.splitext(split)[0] + "_rRNA_Filter_OUT")    
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDSx))
                    break
                PBS_script_out.write(line + "\n")
        JobIDx = subprocess.check_output("ssh gpc01 " + "\"" + "cd " + os.path.dirname(Split_File) + ";" + "qsub" + " " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs") + "\"", shell=True)
        Preprocess_jobs.append(JobIDx.strip("\n"))
        export_filepath = os.path.splitext(split)[0]
end_unpaired_infernal = clock()

start_paired_infernal = clock()
#calls all parts of a split pair 1, 2
for split in sorted(os.listdir(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"))):
    if split.split("_paired_n_contaminants_split_")[0].endswith("1"):
        Split_File1 = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", split.split("_paired_n_contaminants_split_")[0][:-1] + "1" + "_paired_n_contaminants_split_" + split.split("_paired_n_contaminants_split_")[1])
        Split_File2 = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", split.split("_paired_n_contaminants_split_")[0][:-1] + "2" + "_paired_n_contaminants_split_" + split.split("_paired_n_contaminants_split_")[1])

        COMMANDSy = [Python + " " + Filter_rRNA + " " + Split_File1 + " " + os.path.splitext(Split_File1)[0] + "_mRNA.fastq" + " " + os.path.splitext(Split_File1)[0] + "_rRNA.fastq" + " " + Split_File2 + " " + os.path.splitext(Split_File2)[0] + "_mRNA.fastq" + " " + os.path.splitext(Split_File2)[0] + "_rRNA.fastq"]

        with open(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs"), "w") as PBS_script_out:
            for line in PBS_Submit_LowMem.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(split)[0] + "_rRNA_Filter")
                if "ERROR" in line:
                    line = line.replace("ERROR", os.path.splitext(split)[0] + "_rRNA_Filter_ERR")    
                if "OUTPUT" in line:
                    line = line.replace("OUTPUT", os.path.splitext(split)[0] + "_rRNA_Filter_OUT")    
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDSy))
                    break
                PBS_script_out.write(line + "\n")
        JobIDy = subprocess.check_output("ssh gpc01 " + "\"" + "cd " + os.path.dirname(Split_File1) + ";" + "qsub" + " " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs") + "\"", shell=True)
        Preprocess_jobs.append(JobIDy.strip("\n"))

print ":".join(Preprocess_jobs[-10:])
end_paired_infernal = clock()
end_all = clock()

with open(export_filepath + "_rna_split_jobs.txt", 'w+') as profile:
    
    profile.write("rRNA split job\n")
    profile.write( "============================================\n")
    profile.write("total runtime: " + str(end_all - start_all) + "s\n")
    profile.write( "unpaired infernal runtime: " + str(end_unpaired_infernal - start_unpaired_infernal) + "s\n")
    profile.write( "paired infernal runtime: " + str(end_paired_infernal - start_paired_infernal) + "s\n")
    profile.close()
"""