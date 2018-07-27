# This is the high-level mode that calls the infernal program

#!/usr/bin/env python

import sys
import os

import mt_pipe_paths as mpp
import mt_pipe_commands as mpc

def filter_rRNA(split_folder, category, stage_name):
    #convert fastq splits into fasta
    #then put fasta files through infernal.
    dir_list = os.listdir(split_folder)
    job_id_list = []
    job_count = 0
    for item in dir_list:
        seq_file_name = split_folder + item.split('.')[0]
        inner_name = category + "_" + str(job_count) #the name of the individual files
        comm = mpc
        job_id_list.append(comm.create_and_launch(
                            command_list = comm.create_infernal_command(split_folder + inner_name),
                            job_name = stage_name, 
                            inner_name = inner_name#,
                            #run_job = True
                            ) 
                        )
        job_count += 1
        
if __name__ == "__main__":
    #needs the folder containing all of the split fastqs
    #we will convert here
    split_folder = sys.argv[1]
    category = sys.argv[2]
    stage_name = sys.argv[3]
    filter_rRNA(split_folder, category, stage_name)

