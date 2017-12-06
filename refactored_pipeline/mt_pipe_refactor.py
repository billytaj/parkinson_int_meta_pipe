#!/usr/bin/env python
# this pipeline is a mess.
# This program simply makes queue submissions to the cluster, and in doing so, eats a lot of down
import sys
import os
import os.path
import subprocess
import multiprocessing
import mt_pipe_commands as mpcom
import mt_pipe_paths as mpfp

Threads = str(multiprocessing.cpu_count())

run_jobs = False

def main(input_folder, output_folder):
    # constants
    # -----------------------------
    single_mode = 0
    double_mode = 1

    # system vars
    # -------------------------------
    # operating mode:
    # 0: single-ended
    # 1: paired
    # 2: error
    operating_mode = 0
    
    #note: this also needs to support paired and single-ended data
    #input folder is the main location of the dump.
    #
    #for genome in sorted(os.listdir(input_folder)):
    file_list = []
    Network_list = []
    #only seems to look for *1.fastq, and nothing else.  the whole loop is wasting time.  
    raw_sequence_path = input_folder + "/raw_sequences/"
    if not os.path.exists(raw_sequence_path):
        print("No sequences found.  to use pipeline, please place fastq file at:", raw_sequence_path)
        os.makedirs(raw_sequence_path)
        sys.exit()
    else:
        # folder found.  now see if it's single-ended or paired
        # for single-ended, have only 1 file.  if doubled, have both files.  Else, stop
        genome_file_count = len(os.listdir(raw_sequence_path))
        print("number of files:", genome_file_count)
        if(genome_file_count == 1):
            print("OPERATING IN SINGLE-ENDED MODE")
        elif(genome_file_count == 2):
            print("OPERATING IN PAIRED-MODE")
        else:
            print("Too many genome files here.  Get rid of all non-essentials")
            sys.exit()
        
        operating_mode = genome_file_count - 1
        
        # is the file too big?
        # split it.
        
        #init a command object, and start making commands
        #sys.exit()
        
        if(operating_mode == double_mode):
            preprocess_label = "preprocess"
            raw_pair_0_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[0]
            raw_pair_1_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[1]
            comm = mpcom.mt_pipe_commands(Quality_score = 33, Thread_count = 16, raw_sequence_path_0 = raw_pair_0_path, raw_sequence_path_1 = raw_pair_1_path)
            #preprocess_job_id = comm.create_pbs_and_launch("preprocess", comm.create_pre_double_command(preprocess_label), run_job = True)
            #rRNA_job_id = comm.create_pbs_and_launch("rRNA_filter", comm.create_rRNA_filter_command("rRNA_filter", 5, "preprocess"), dependency_list = preprocess_job_id, run_job = True)
            rRNA_job_id = comm.create_pbs_and_launch("rRNA_filter", comm.create_rRNA_filter_command("rRNA_filter", 5, "preprocess"), run_job = True)
            
        elif(operating_mode == single_mode):
            print("not ready")

        """
            # Preprocessing
            

            create_pbs_job("Preprocess", Input_FName, COMMANDS_Pre)        
            if(run_jobs):
                JobID_Pre = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Preprocess.pbs"])

            
            create_pbs_job("rRNA_Submit", Input_FName, COMMANDS_rRNA)
            if(run_jobs):
                JobID_rRNA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_rRNA_Submit.pbs", "-W", "depend=afterok:" + JobID_Pre.strip("\n")])


            create_pbs_job("Combine", Input_FName, COMMANDS_Combine)
            if(run_jobs):
                JobID_Combine = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine.pbs", "-W", "depend=afterok:" + JobID_rRNA.strip("\n")])

                subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine.strip("\n").split(".")[0], JobID_rRNA.strip("\n")])

            #BIGG Database, AGORA Nature paper, Additional functionality
            # Transcript Assembly
            Contigs = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_SpadesOut", "contigs.fasta")

            create_pbs_job("Assemble", Input_FName, COMMANDS_Assemble)        
            if(run_jobs):
                JobID_Assemble = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Assemble.pbs", "-W", "depend=afterok:" + JobID_Combine.strip("\n")])

            # Protein Annotation BWA
            create_pbs_job("Annotate_BWA", Input_FName, COMMANDS_Annotate_BWA, "med")        
            if(run_jobs):
                JobID_Annotate_BWA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BWA.pbs", "-W", "depend=afterok:" + JobID_Assemble.strip("\n")])

            # Protein Annotation BLAT 1
            create_pbs_job("Annotate_BLAT1", Input_FName, COMMANDS_Annotate_BLAT1, "med")        
            if(run_jobs):
                JobID_Annotate_BLAT1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 2
            create_pbs_job("Annotate_BLAT2", Input_FName, COMMANDS_Annotate_BLAT2, "med")
            if(run_jobs):
                JobID_Annotate_BLAT2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 3
            create_pbs_job("Annotate_BLAT3", Input_FName, COMMANDS_Annotate_BLAT3, "med")       
            if(run_jobs):
                JobID_Annotate_BLAT3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 4
            create_pbs_job("Annotate_BLAT4", Input_FName, COMMANDS_Annotate_BLAT4, "med")
            if(run_jobs):
                JobID_Annotate_BLAT4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT Postprocessing
            create_pbs_job("Annotate_BLAT_Postprocessing", Input_FName, COMMANDS_Annotate_BLAT_Post)
            if(run_jobs):
                JobID_Annotate_BLAT_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT_Postprocessing.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT1.strip("\n") + ":" + JobID_Annotate_BLAT2.strip("\n") + ":" + JobID_Annotate_BLAT3.strip("\n") + ":" + JobID_Annotate_BLAT4.strip("\n")])


            # Protein Annotation Diamond 1
            create_pbs_job("Annotate_DMD", Input_FName, COMMANDS_Annotate_Diamond1)        
            if(run_jobs):
                JobID_Annotate_Diamond1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 2
            create_pbs_job("Annotate_DMD2", Input_FName, COMMANDS_Annotate_Diamond2)  
            if(run_jobs):
                JobID_Annotate_Diamond2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 3
            create_pbs_job("Annotate_DMD3", Input_FName, COMMANDS_Annotate_Diamond3)
            if(run_jobs):
                JobID_Annotate_Diamond3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 4
          
            create_pbs_job("Annotate_DMD4", Input_FName, COMMANDS_Annotate_Diamond4)                
            if(run_jobs):
                JobID_Annotate_Diamond4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond Postprocess
            create_pbs_job("Annotate_DMD_Postprocessing", Input_FName, COMMANDS_Annotate_Diamond_Post)
            if(run_jobs):
                JobID_Annotate_Diamond_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD_Postprocess.pbs", "-W", "depend=afterany:" + JobID_Annotate_Diamond1.strip("\n") + ":" + JobID_Annotate_Diamond2.strip("\n") + ":" + JobID_Annotate_Diamond3.strip("\n") + ":" + JobID_Annotate_Diamond4.strip("\n")])
            
            # Classify Reads
            create_pbs_job("Classify", Input_FName, COMMANDS_Classify, "high")
            if(run_jobs):
                JobID_Classify = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Classify.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])
            
            # Prepare EC annotation files
            EC_Split = os.path.join(Input_Filepath + "_EC_Annotation", "Split")
            EC_Output = os.path.join(Input_Filepath + "_EC_Annotation", "Output")
            
            create_pbs_job("EC_Preprocess", Input_FName, COMMANDS_EC_Preprocess)
            if(run_jobs):
                JobID_EC_Preprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Preprocess.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])

            #EC detection
            create_pbs_job("Detect", Input_FName, COMMANDS_Detect)        
            if(run_jobs):
                JobID_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Detect.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            #Combine the EC files
            create_pbs_job("Combine_Detect", Input_FName, COMMANDS_Combine_Detect)
            if(run_jobs):
                JobID_Combine_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine_Detect.pbs", "-W", "depend=afterok:" + JobID_Detect.strip("\n")])
            
                subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine_Detect.strip("\n").split(".")[0], JobID_Detect.strip("\n")])

            
            #PRIAM stage
            create_pbs_job("PRIAM", Input_FName, COMMANDS_PRIAM)        
            if(run_jobs):
                JobID_PRIAM = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_PRIAM.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            # EC Diamond
            create_pbs_job("EC_Diamond", Input_FName, COMMANDS_EC_Diamond)        
            if(run_jobs):
                JobID_EC_Diamond = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Diamond.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            # EC Annotation Compile
            

            with open(os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_EC_Postprocess")
                    if "ERROR" in line:
                        line = line.replace("ERROR", os.path.splitext(Input_FName)[0] + "_EC_Postprocess_ERR")
                    if "OUTPUT" in line:
                        line = line.replace("OUTPUT", os.path.splitext(Input_FName)[0] + "_EC_Postprocess_OUT")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_EC_Postprocess))
                        break
                    PBS_script_out.write(line + "\n")
            if(run_jobs):
                JobID_EC_Postprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "-W", "depend=afterok:" + JobID_Combine_Detect.strip("\n") + ":" + JobID_PRIAM.strip("\n") + ":" + JobID_EC_Diamond.strip("\n")])

            # Network Generation
            

            with open(os.path.splitext(Input_FName)[0] + "_Network.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Network")
                    if "ERROR" in line:
                        line = line.replace("ERROR", os.path.splitext(Input_FName)[0] + "_Network_ERR")
                    if "OUTPUT" in line:
                        line = line.replace("OUTPUT", os.path.splitext(Input_FName)[0] + "_Network_OUT")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Network))
                        break
                    PBS_script_out.write(line + "\n")
            
            if(run_jobs):
                JobID_Network = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Network.pbs", "-W", "depend=afterok:" + JobID_EC_Postprocess.strip("\n") + ":" + JobID_Classify.strip("\n")])
            Network_list.append(JobID_Network.strip("\n"))

    if len(file_list) > 1:
        os.chdir(output_folder)
        Input_Filepath = os.path.join(output_folder, os.path.splitext(os.path.basename(file_list[0]))[0])[:-11]
        
        with open(os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "w") as PBS_script_out:
            for line in PBS_Submit_LowMem.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join")
                if "ERROR" in line:
                    line = line.replace("ERROR", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join_ERR")
                if "OUTPUT" in line:
                    line = line.replace("OUTPUT", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join_OUT")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDS_Join))
                    break
                PBS_script_out.write(line + "\n")
        
        if(run_jobs):
            JobID_Join = subprocess.check_output(["qsub", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "-W", "depend=afterok:" + ":".join(Network_list)])
"""                

if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    os.chdir(output_folder)
    main(input_folder, output_folder)