#!/usr/bin/env python
# this pipeline is a mess.
# This program simply makes queue submissions to the cluster, and in doing so, eats a lot of down
import sys
import os
import os.path
import subprocess as sp
import multiprocessing as mp
import mt_pipe_commands as mpcom
import mt_pipe_paths as mpfp
import time

run_jobs = False
def make_folder(folder_path):
    if not(os.path.exists(folder_path)):
        os.makedirs(folder_path)

class qsub_sync:
    def __init__(self, system_mode):
        self.system_mode = system_mode

    #handles where to auto-resume
    def check_where_resume(self, job_label = None, full_path = None):
        if(job_label):
            job_path = job_label + "/data/final_results"
            print("looking at:", job_path)
            if(os.path.exists(job_path)):
                file_list = os.listdir(job_path)
                if(len(file_list) > 0):
                    print("bypassing!")
                    return True
                else:
                    print("running")
                    return False
                    
            else:
                print("doesn't exist: running")
                return False
        else:
            job_path = full_path
            print("looking at:", job_path)
            if(os.path.exists(job_path)):
                file_list = os.listdir(job_path)
                if(len(file_list) > 0):
                    print("bypassing!")
                    return True
                else:
                    print("running")
                    return False
                    
            else:
                print("doesn't exist: running")
                return False
    #handles job syncing
    # op mode:  which job category -> completed, running, in queue, blocked
    # check_mode: how to check our job IDs vs the job category list 
    # -> any: true on existence of 1
    # -> all: true on existence of all 
    
    def check_job_finished(self, job_id, op_mode, check_mode="all"):
        jobs_list = []
        try:
            if(op_mode == "completed"):
                jobs_list = sp.check_output(["showq", "-u billyc59", "-c"]).decode('ascii') #completed
            elif(op_mode == "queue"):
                jobs_list = sp.check_output(["showq", "-u billyc59", "-i"]).decode('ascii') #in the queue
            elif(op_mode == "run"):
                jobs_list = sp.check_output(["showq", "-u billyc59", "-r"]).decode('ascii') #running
            elif(op_mode == "blocked"):
                jobs_list = sp.check_output(["showq", "-u billyc59", "-b"]).decode('ascii') #holding
            else:
                print("bad op mode to check_job_finished.  killing")
                sys.exit()
               
            #print(op_mode)
            #print("Jobs list from qstat:", jobs_list)
            if (isinstance(job_id, int)):
                #single
                if(str(job_id) in jobs_list):
                    return 1 #match
                else:
                    return 2 #no match
                
            elif(isinstance(job_id, list)):
                # multi.  
                if(check_mode == "all"):
                    #all: return true only if all are in the list
                    match_count = 0
                    for item in job_id:
                        if(str(item) in jobs_list):
                            match_count += 1
                    if(match_count != len(job_id)):
                        return 2 #no match
                    else:
                        return 1 #match
                        
                elif(check_mode == "any"):
                    #any: true on any hit of job id matching the list contents
                    #print("-----------------------------------")
                    for item in job_id:
                        #print("looking for:", str(item))
                        if(str(item) in jobs_list):
                            #print("key found:", str(item))
                            return 1 #found
                        #else:
                            #print("can't find key in list")
                    return 2 #not found
                
            elif(job_id is None):
                print("bad arg to check job finished")
                return 0 #error
            else:
                print("This isn't supposed to happen")
                sys.exit()
        except Exception as e:
            print("exception happened with syncer:", e)
            print("we're ending everything")
            sys.exit()
            
    def wait_for_sync(self, timeout, job_id, label, message = None):
        #does the actual wait for the qsub job to finish
        #no longer needed, with multiprocessing -> join
        if((self.system_mode == "docker") or (self.system_mode == "singularity")):
            print("running in docker.  bypassing sync")
        else:
            if(job_id is None) or (not job_id):
                print("Blank job id.  killing")
                sys.exit()
            
            #print("job id list:")
            #print(job_id)
            time.sleep(60)
            print(message)
            b_lock = True
            job_status = 0
            message_instance_timeout = 0
            lockout_count = int(timeout)
            
            retry_counter = 10    
            while(b_lock):
                job_status = self.check_job_finished(job_id, "completed", "all")
                if(job_status == 1):
                    b_lock = False
                else:
                    
                    #only start the countdown if the job is running, not while it's waiting in the queue.  
                    #job could wait a long time in the queue if someone's clobbering the cluster
                    job_status = self.check_job_finished(job_id, "run", "any")
                    if(job_status == 1):
                        lockout_count -= 1
                        #print("job still running")
                    elif(job_status == 0):
                        #print("something happened.  we're dying")
                        sys.exit()
                    elif(job_status == 2):
                        #print("not running")
                        #none running.  maybe in queue?
                        job_status = self.check_job_finished(job_id, "queue", "any")
                        if(job_status == 2):
                            #print("not in queue")
                            #nothing in queue.  maybe blocked?
                            job_status = self.check_job_finished(job_id, "blocked", "any")
                            if(job_status == 2):
                                #print("not blocked")
                                #not in any queue
                                #print("sync over")
                                if(retry_counter == 0):
                                    b_lock = False
                                else:
                                    retry_counter -= 1
                        #    elif(job_status == 1):
                        #        print("job blocked")
                        #elif(job_status == 1):
                        #    print("job in queue")
                    
                    else:
                        print("this can't happen in sync.  killing")
                        sys.exit()
                
                if(lockout_count <= 0):
                    print(label, "took too long.  shutting down pipeline.  QSUB Jobs may still be running")
                    print("Use qselect -u <username> | xargs qdel  to remove jobs from qsub")
                    sys.exit()
                
                time.sleep(5)
                
            
def main(input_folder, output_folder, system_op):
    
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
    sync_obj = qsub_sync(system_op)
    start_time = time.time()
    #note: this also needs to support paired and single-ended data
    #input folder is the main location of the dump.
    #
    #for genome in sorted(os.listdir(input_folder)):
    file_list = []
    Network_list = []
    
    
    mp_store = [] #stores the multiprocessing processes
    
    #only seems to look for *1.fastq, and nothing else.  the whole loop is wasting time.  
    raw_sequence_path = input_folder #+ "/raw_sequences/"
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
        #operating mode denotes whether it's single or paired reads
        operating_mode = genome_file_count - 1 
        
        # is the file too big?
        # split it.
        
        #init a command object, and start making commands
        #sys.exit()
        
        
        
        
        
        if(operating_mode == double_mode):
            preprocess_label = "preprocess"
            rRNA_filter_label = "rRNA_filter"
            repop_job_label = "duplicate_repopulation"
            assemble_contigs_label = "assemble_contigs"
            gene_annotation_BWA_label = "gene_annotation_BWA"
            gene_annotation_BLAT_label = "gene_annotation_BLAT"
            GA_BLAT_PP_label = "BLAT_postprocess"
            gene_annotation_DIAMOND_label = "gene_annotation_DIAMOND"
            
            rRNA_filter_orphans_fastq_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_fastq/"
            rRNA_filter_pair_1_fastq_folder = os.getcwd()  + "/rRNA_filter/data/pair_1/pair_1_fastq/"
            rRNA_filter_pair_2_fastq_folder = os.getcwd()  + "/rRNA_filter/data/pair_2/pair_2_fastq/"
                        
            raw_pair_0_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[0]
            raw_pair_1_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[1]
            comm = mpcom.mt_pipe_commands(Quality_score = 33, Thread_count = 16, system_mode = system_op, raw_sequence_path_0 = raw_pair_0_path, raw_sequence_path_1 = raw_pair_1_path) #start obj
            
            #need to erase the job ids.  they're no longer needed
            if(not sync_obj.check_where_resume(output_folder + preprocess_label)):
                
                process = mp.Process(
                    target = comm.create_pbs_and_launch, 
                    args = (
                        preprocess_label, 
                        comm.create_pre_double_command(preprocess_label), 
                        True
                    )
                )
                process.start()
                process.join()
                #preprocess_job_id = None
            
            #else:
            #    preprocess_job_id = None
            
            #rRNA_filter_job_id = []
            
            if(not sync_obj.check_where_resume(output_folder +  rRNA_filter_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch, 
                    args = (
                        rRNA_filter_label, 
                        comm.create_rRNA_filter_prep_command(
                        rRNA_filter_label, 5, preprocess_label), 
                        #dependency_list = preprocess_job_id, 
                        True
                    )
                )
                process.start()
                process.join()
                orphans_mRNA_path = output_folder + rRNA_filter_label + "/data/orphans/orphans_mRNA"
                if(not sync_obj.check_where_resume(None, orphans_mRNA_path)):
                    for item in os.listdir(rRNA_filter_orphans_fastq_folder):
                        file_root_name = item.split('.')[0]
                        inner_name = file_root_name + "_infernal"
                        process = mp.Process(
                            target = comm.create_pbs_and_launch, 
                            args = (
                                "rRNA_filter", 
                                comm.create_rRNA_filter_command("rRNA_filter", "orphans", file_root_name), 
                                
                                #dependency_list = rRNA_filter_job_id[0],
                                True,
                                inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process)
                #pausing here for now, because we want to stage the infernal calls
                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list
                
                
                for item in os.listdir(rRNA_filter_pair_1_fastq_folder):
                    file_root_name = item.split('.')[0]
                    inner_name = file_root_name + "_infernal"
                    process = mp.Process(
                        target = comm.create_pbs_and_launch,
                        args = (
                            "rRNA_filter", 
                            comm.create_rRNA_filter_command("rRNA_filter", "pair_1", file_root_name), 
                            #dependency_list = rRNA_filter_job_id[0],
                            True,
                            inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
                    
                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list
                    
                for item in os.listdir(rRNA_filter_pair_2_fastq_folder):
                    file_root_name = item.split('.')[0]
                    inner_name = file_root_name + "_infernal"
                    process = mp.Process(
                        target = comm.create_pbs_and_launch,
                        args = (
                            "rRNA_filter", 
                            comm.create_rRNA_filter_command("rRNA_filter", "pair_2", file_root_name), 
                            #dependency_list = rRNA_filter_job_id[0],
                            True,
                            inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list
                
                #wait for infernal to finish running
                #time.sleep(5)
                #sync_obj.wait_for_sync(800, rRNA_filter_job_id, rRNA_filter_label, "waiting for Infernal")
                #print("rRNA ID list:", rRNA_filter_job_id)
                inner_name = "rRNA_filter_post"
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                        rRNA_filter_label, 
                        comm.create_rRNA_filter_post_command(rRNA_filter_label), 
                        True,
                        inner_name
                    )
                )
                process.start()
                process.join()
                
            #print("ending prematurely")
            #sys.exit()
            
            #-------------------------------------------------------------
            #Next, we have duplicate repopulation
            if(not sync_obj.check_where_resume(output_folder +  repop_job_label)):
            
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    repop_job_label, 
                    comm.create_repop_command(repop_job_label, preprocess_label, rRNA_filter_label), 
                    True
                    )
                )
                process.start()
                process.join()
            
            #----------------------------------------
            # assemble contigs
            if(not sync_obj.check_where_resume(output_folder + assemble_contigs_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    assemble_contigs_label, 
                    comm.create_assemble_contigs_command(assemble_contigs_label, repop_job_label),
                    True
                    )
                )
                process.start()
                process.join()
            #else:
            #    assemble_contigs_id = None
            
            
            #----------------------------------------------
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_BWA_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_BWA_label,
                    comm.create_BWA_annotate_command(gene_annotation_BWA_label, assemble_contigs_label),
                    True
                    )
                )
                process.start()
                process.join()
            #else:
            #    gene_annotation_BWA_id = None
                
                
            #------------------------------------------------
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_BLAT_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_BLAT_label, 
                    comm.create_BLAT_annotate_command( gene_annotation_BLAT_label, gene_annotation_BWA_label),
                    #dependency_list = gene_annotation_BWA_id,
                    True
                    )
                )
                process.start()
                process.join()
                
            
            #------------------------------------------------
            if(not sync_obj.check_where_resume(output_folder + GA_BLAT_PP_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    GA_BLAT_PP_label, 
                    comm.create_BLAT_pp_command(GA_BLAT_PP_label, gene_annotation_BWA_label, gene_annotation_BLAT_label),
                    True
                    )
                )
                process.start()
                process.join()
                
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_DIAMOND_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_DIAMOND_label,
                    comm.create_DIAMOND_annotate_command(gene_annotation_DIAMOND_label, GA_BLAT_PP_label),
                    True
                    )
                )
                process.start()
                process.join()
                
            end_time = time.time()
            print("Total runtime:", end_time - start_time)
            
            
            
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
    if(len(sys.argv) < 3):
        print("no args provided.  try again:  arg(1) input folder, arg(2) output folder, arg(3) docker or scinet")
        sys.exit()
    else:    
        input_folder = sys.argv[1]
        output_folder = sys.argv[2]
        system_op = sys.argv[3] #scinet or docker
        #scinet_user_name = sys.argv[3]
        if not(os.path.exists(output_folder)):
            print("output file doesn't exist.  now making one")
            os.makedirs(output_folder)
        os.chdir(output_folder)
        main(input_folder, output_folder, system_op)