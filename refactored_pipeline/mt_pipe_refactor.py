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

class sync_control:
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


def main(input_folder, output_folder, system_op, user_mode):

    # constants
    # -----------------------------
    single_mode = 1
    double_mode = 2

    # system vars
    # -------------------------------
    # operating mode:
    # 0: single-ended
    # 1: paired
    # 2: error
    operating_mode = 0
    sync_obj = sync_control(system_op)
    thread_count = mp.cpu_count() #should not be hard-coded
    
    #note: this also needs to support paired and single-ended data
    #input folder is the main location of the dump.
    #
    file_list = []
    Network_list = []


    mp_store = [] #stores the multiprocessing processes
    
    #--------------------------------------------------
    #profiling vars
    start_time              = 0
    end_time                = 0
    preprocess_start        = 0
    preprocess_end          = 0
    rRNA_filter_start       = 0
    rRNA_filter_end         = 0
    repop_start             = 0
    repop_end               = 0
    assemble_contigs_start  = 0
    assemble_contigs_end    = 0
    GA_BWA_start            = 0
    GA_BWA_end              = 0
    GA_BLAT_start           = 0
    GA_BLAT_end             = 0
    GA_DIAMOND_start        = 0
    GA_DIAMOND_end          = 0
    TA_start                = 0
    TA_end                  = 0
    EC_start                = 0
    EC_end                  = 0
    EC_DETECT_start         = 0
    EC_DETECT_end           = 0
    EC_PRIAM_DIAMOND_start  = 0
    EC_PRIAM_DIAMOND_end    = 0
    Cytoscape_start         = 0
    Cytoscape_end           = 0
    
    start_time = time.time()
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
        operating_mode = genome_file_count
        if(genome_file_count == 0):
            print("There's no genomes to analyze.  we're shutting down")
            sys.exit()
        elif(genome_file_count == 1):
            print("OPERATING IN SINGLE-ENDED MODE")
        elif(genome_file_count == 2):
            print("OPERATING IN PAIRED-MODE")
        else:
            print("Too many genome files here.  Get rid of all non-essentials")
            sys.exit()
        #operating mode denotes whether it's single or paired reads


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
            GL_BLAT_cat_label = "BLAT_cat_results"
            GA_BLAT_PP_label = "BLAT_postprocess"
            gene_annotation_DIAMOND_label = "gene_annotation_DIAMOND"
            taxon_annotation_label = "taxonomic_annotation"
            ec_annotation_label = "enzyme_annotation"
            network_label = "RPKM_network"

            rRNA_filter_orphans_fastq_folder = os.getcwd() + "/rRNA_filter/data/orphans/orphans_fastq/"
            rRNA_filter_pair_1_fastq_folder = os.getcwd()  + "/rRNA_filter/data/pair_1/pair_1_fastq/"
            rRNA_filter_pair_2_fastq_folder = os.getcwd()  + "/rRNA_filter/data/pair_2/pair_2_fastq/"

            raw_pair_0_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[0]
            raw_pair_1_path = raw_sequence_path + sorted(os.listdir(raw_sequence_path))[1]
            quality_encoding = 33 #This should not be a constant, we need some process to determine the quality encoding ex. vsearch --fastq_chars
            comm = mpcom.mt_pipe_commands(Quality_score = quality_encoding, Thread_count = thread_count, system_mode = system_op, user_mode = user_mode, raw_sequence_path_0 = raw_pair_0_path, raw_sequence_path_1 = raw_pair_1_path) #start obj

            #Creates our command object, for creating shellscripts.
            comm = mpcom.mt_pipe_commands(
                                            Quality_score = quality_encoding,                     #leftover from an argument into one of the tools
                                            Thread_count = thread_count,                      #leftover, same reason.
                                            system_mode = system_op,                #docker, or singularity
                                            raw_sequence_path_0 = raw_pair_0_path,  #
                                            raw_sequence_path_1 = raw_pair_1_path
                                        ) #start obj

            #This is the format we use to launch each stage of the pipeline.
            #We start a multiprocess that starts a subprocess.
            #The subprocess is created from the comm object

            #The preprocess stage
            preprocess_start = time.time()
            if(not sync_obj.check_where_resume(output_folder + preprocess_label)):

                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                        preprocess_label,
                        comm.create_pre_double_command(preprocess_label),
                        True
                    )
                )
                process.start() #start the multiprocess
                process.join()  #wait for it to end
            preprocess_end = time.time()
            
            #rRNA removal stage
            rRNA_filter_start = time.time()
            if(not sync_obj.check_where_resume(output_folder +  rRNA_filter_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                        rRNA_filter_label,
                        comm.create_rRNA_filter_prep_command(
                        rRNA_filter_label, int(mp.cpu_count()/2), preprocess_label),
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
                                True,
                                inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process) # pack all the processes into a list
                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list

                pair_1_mRNA_path = output_folder + rRNA_filter_label + "/data/pair_1/pair_1_mRNA"
                if(not sync_obj.check_where_resume(None, pair_1_mRNA_path)):
                    for item in os.listdir(rRNA_filter_pair_1_fastq_folder):
                        file_root_name = item.split('.')[0]
                        inner_name = file_root_name + "_infernal"
                        process = mp.Process(
                            target = comm.create_pbs_and_launch,
                            args = (
                                "rRNA_filter",
                                comm.create_rRNA_filter_command("rRNA_filter", "pair_1", file_root_name),
                                True,
                                inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process)
                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list

                pair_2_mRNA_path = output_folder + rRNA_filter_label + "/data/pair_2/pair_2_mRNA"
                if(not sync_obj.check_where_resume(None, pair_2_mRNA_path)):
                    for item in os.listdir(rRNA_filter_pair_2_fastq_folder):
                        file_root_name = item.split('.')[0]
                        inner_name = file_root_name + "_infernal"
                        process = mp.Process(
                            target = comm.create_pbs_and_launch,
                            args = (
                                "rRNA_filter",
                                comm.create_rRNA_filter_command("rRNA_filter", "pair_2", file_root_name),
                                True,
                                inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process)

                for item in mp_store:
                    item.join() # wait for things to finish
                mp_store[:] = [] #clear the list

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

            rRNA_filter_end = time.time()
            #-------------------------------------------------------------
            # Duplicate repopulation
            repop_start = time.time()
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
            repop_end = time.time()
            #----------------------------------------
            # Assemble contigs
            assemble_contigs_start = time.time()
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
            assemble_contigs_end = time.time()

            #----------------------------------------------
            # BWA gene annotation
            GA_BWA_start = time.time()
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_BWA_label)):

                names = ["contigs", "orphans", "pair_1", "pair_2"]
                mp_store[:] = []
                for item in names:
                    inner_name = "BWA_" + item
                    process = mp.Process(
                        target = comm.create_pbs_and_launch,
                        args = (
                        gene_annotation_BWA_label,
                        comm.create_BWA_annotate_command(gene_annotation_BWA_label, assemble_contigs_label, item),
                        True,
                        inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
                
                for item in mp_store:
                    item.join()
                mp_store[:] = [] #clear the list

            #else:
            #    gene_annotation_BWA_id = None
                inner_name = "BWA_pp"
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_BWA_label,
                    comm.create_BWA_pp_command(gene_annotation_BWA_label, assemble_contigs_label),
                    True,
                    inner_name
                    )
                )
                process.start()
                process.join()
            GA_BWA_end = time.time()
            #------------------------------------------------
            # BLAT gene annotation
            GA_BLAT_start = time.time()
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_BLAT_label)):

                split = 5#mp.cpu_count() #split based on the way microbial_cds_db was split.  this must also change
                names = ["contigs", "orphans", "pair_1", "pair_2"]
                for item in names:
                    for i in range (1, split+1):
                        inner_name = "BLAT_"+item + "_"+str(i)
                        process = mp.Process(
                            target = comm.create_pbs_and_launch,
                            args = (
                            gene_annotation_BLAT_label, 
                            comm.create_BLAT_annotate_command(gene_annotation_BLAT_label, gene_annotation_BWA_label, item, i),
                            #dependency_list = gene_annotation_BWA_id,
                            True,
                            inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process)
                for item in mp_store:        
                    item.join()
                mp_store[:] = [] #clear the list

                for item in names:
                    inner_name = item + "_cat"
                    process = mp.Process(
                        target = comm.create_pbs_and_launch,
                        args = (
                        gene_annotation_BLAT_label, 
                        comm.create_BLAT_cat_command(gene_annotation_BLAT_label, item, split),
                        True,
                        inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
                for item in mp_store:
                    item.join()
                mp_store[:] = []

                inner_name = "BLAT_pp"
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_BLAT_label, 
                    comm.create_BLAT_pp_command(gene_annotation_BLAT_label, gene_annotation_BWA_label),
                    True,
                    inner_name
                    )
                )
                process.start()
                process.join()
            GA_BLAT_end = time.time()

            #------------------------------------------------------
            # Diamond gene annotation
            GA_DIAMOND_start = time.time()
            if(not sync_obj.check_where_resume(output_folder + gene_annotation_DIAMOND_label)):

                names = ["contigs", "orphans", "pair_1", "pair_2"]
                for item in names:
                    inner_name = item + "_run_diamond"
                    process = mp.Process(
                        target = comm.create_pbs_and_launch,
                        args = (
                        gene_annotation_DIAMOND_label,
                        comm.create_DIAMOND_annotate_command(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label, item),
                        True,
                        inner_name
                        )
                    )
                    process.start()
                    mp_store.append(process)
                for item in mp_store:
                    item.join()
                mp_store[:] = []

                inner_name = "diamond_pp"
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    gene_annotation_DIAMOND_label,
                    comm.create_DIAMOND_pp_command(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label),
                    True,
                    inner_name
                    )
                )
                process.start()
                process.join()
            GA_DIAMOND_end = time.time()
            # ------------------------------------------------------
            # Taxonomic annotation
            TA_start = time.time()
            if(not sync_obj.check_where_resume(output_folder + taxon_annotation_label)):
                process = mp.Process(
                    target = comm.create_pbs_and_launch,
                    args = (
                    taxon_annotation_label,
                    comm.create_taxonomic_annotation_command(taxon_annotation_label, assemble_contigs_label, gene_annotation_DIAMOND_label),
                    True
                    )
                )
                process.start()
                process.join()
            TA_end = time.time()
            
            # ------------------------------------------------------
            # EC annotation
            EC_start = time.time()
            EC_DETECT_start = time.time()
            if (not sync_obj.check_where_resume(output_folder + ec_annotation_label)):
                # Preparing folders for DETECT
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        taxon_annotation_label,
                        comm.create_EC_DETECT_prep(ec_annotation_label, gene_annotation_DIAMOND_label, int(mp.cpu_count()/2)),
                        True
                    )
                )
                process.start()
                process.join()
 
                # Running DETECT on split protien files 
                proteins_path = output_folder + ec_annotation_label + "/data/0_proteins/"
                if (not sync_obj.check_where_resume(None, proteins_path)):
                    for item in os.listdir(proteins_path):
                        file_root_name = os.path.splitext(item)[0]
                        inner_name = file_root_name + "_detect"
                        process = mp.Process(
                            target=comm.create_pbs_and_launch,
                            args=(
                                ec_annotation_label,
                                comm.create_EC_DETECT_command(ec_annotation_label, file_root_name),
                                True,
                                inner_name
                            )
                        )
                        process.start()
                        mp_store.append(process)  # pack all the processes into a list
 
                for item in mp_store:
                    item.join()  # wait for things to finish
                mp_store[:] = []  # clear the list
            EC_DETECT_end = time.time()
            DETECT_path = output_folder + ec_annotation_label + "/data/1_detect/"

            # Running Priam and Diamond
            EC_PRIAM_DIAMOND_start = time.time()
            if (not sync_obj.check_where_resume(output_folder + DETECT_path)):
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        ec_annotation_label,
                        comm.create_EC_PRIAM_DIAMOND_command(ec_annotation_label, assemble_contigs_label, gene_annotation_DIAMOND_label),
                        True
                    )
                )
                process.start()
                process.join()
 
                inner_name = "ea_post"
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        ec_annotation_label,
                        comm.create_EC_postprocess_command(ec_annotation_label),
                        True,
                        inner_name
                    )
                )
                process.start()
                process.join()
            EC_PRIAM_DIAMOND_end = time.time()
            EC_end = time.time()
            # ------------------------------------------------------
            # RPKM Table and Cytoscape Network
            Cytoscape_start = time.time()
            if (not sync_obj.check_where_resume(output_folder + network_label)):
                process = mp.Process(
                    target=comm.create_pbs_and_launch,
                    args=(
                        network_label,
                        comm.create_Network_generation_command(network_label, gene_annotation_DIAMOND_label, taxon_annotation_label, ec_annotation_label),
                        True
                    )
                )
                process.start()
                process.join()
            
            Cytoscape_end = time.time()
            end_time = time.time()
            print("Total runtime:", end_time - start_time, "s", "start:", start_time, "end:", end_time)
            print("preprocess:", preprocess_end - preprocess_start, "s", "start:", preprocess_start, "end:", preprocess_end)
            print("rRNA filter:", rRNA_filter_end - rRNA_filter_start, "s", "start:", rRNA_filter_start, "end:", rRNA_filter_end)
            print("repop:", repop_end - repop_start, "s", "start:", repop_start, "end:", repop_end)
            print("assemble contigs:", assemble_contigs_end - assemble_contigs_start, "s", "start:", assemble_contigs_start, "end:", assemble_contigs_end)
            print("GA BWA:", GA_BWA_end - GA_BWA_start, "s", "start:", GA_BWA_start, "end:", GA_BWA_end)
            print("GA BLAT:", GA_BLAT_end - GA_BLAT_start, "s", "start:", GA_BLAT_start, "end:", GA_BLAT_end)
            print("GA DIAMOND:", GA_DIAMOND_end - GA_DIAMOND_start, "s", "start:", GA_DIAMOND_start, "end:", GA_DIAMOND_end)
            print("TA:", TA_end - TA_start, "s", "start:", TA_start, "end:", TA_end)
            print("EC:", EC_end - EC_start, "s", "start:", EC_start, "end:", EC_end)
            print("EC DETECT:", EC_DETECT_end - EC_DETECT_start, "s", "start:", EC_DETECT_start, "end:", EC_DETECT_end)
            print("EC PRIAM + DIAMOND:", EC_PRIAM_DIAMOND_end - EC_PRIAM_DIAMOND_end, "s", "start:", EC_PRIAM_DIAMOND_start, "end:", EC_PRIAM_DIAMOND_end)
            print("Cytoscape:", Cytoscape_end - Cytoscape_start, "s", "start:", Cytoscape_start, "end:", Cytoscape_end)

        elif(operating_mode == single_mode):
            print("not ready")

        

if __name__ == "__main__":
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    if(len(sys.argv) < 4):
        print("no args provided.  try again:")
        print("arg(1) input folder")
        print("arg(2) output folder")
        print("arg(3) docker or singularity")
        print("arg(4): billy or bj: for pipeline path location modes")
        sys.exit()
    else:    
        input_folder = sys.argv[1]
        if(not input_folder.endswith("/")):
            input_folder += "/"
        output_folder = sys.argv[2]
        if(not output_folder.endswith("/")):
            output_folder += "/"
            
        system_op = sys.argv[3] #scinet or docker
        #scinet_user_name = sys.argv[3]
        if not(os.path.exists(output_folder)):
            print("output file doesn't exist.  now making one")
            os.makedirs(output_folder)
        os.chdir(output_folder)
        user_mode = sys.argv[4]
        main(input_folder, output_folder, system_op, user_mode)