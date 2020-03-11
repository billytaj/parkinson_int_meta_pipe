#!/usr/bin/env python
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu

def mem_checker(threshold):
    #threshold is a percentage
    mem = psu.virtual_memory()
    available_mem = mem.available
    total_mem = mem.total
    
    available_pct = 100 * available_mem / total_mem
    
    if(available_pct <= threshold):
        return False
    else:
        return True
        
        
def make_folder(folder_path):
    if not (os.path.exists(folder_path)):
        os.makedirs(folder_path)

def delete_folder(folder_path):
    if (os.path.exists(os.path.join(folder_path, "data"))):
        print("deleting", os.path.join(folder_path, "data"))
        shutil.rmtree(os.path.join(folder_path, "data"))
        
def compress_folder(folder_path):
    zip_loc = os.path.join(folder_path, "data")
    z = zipfile.ZipFile(folder_path + "_data.zip", "a", zipfile.ZIP_DEFLATED)
    print("compressing interim files:", folder_path)
    for root, dirs, files in os.walk(zip_loc):
        #print("root:", root)
        #print("dirs:", dirs)
        #print("files:", files)
        #print("===============================")
        for file in files:
            z.write(os.path.join(root, file))
    z.close()
        
def write_to_bypass_log(folder_path, message):
    bypass_log_path = os.path.join(folder_path, "bypass_log.txt")
    with open(bypass_log_path, "a") as bypass_log:
        bypass_log.write("\n")
        new_message = message + "\n"
        bypass_log.write(new_message)
        


def check_bypass_log(folder_path, message):
    bypass_keys_list = list()
    bypass_log_path = os.path.join(folder_path, "bypass_log.txt")
    if(os.path.exists(bypass_log_path)):
        with open(bypass_log_path, "r") as bypass_log:
            for line in bypass_log:
                bypass_key = line.strip("\n")
                bypass_keys_list.append(bypass_key)
        if(message in bypass_keys_list):
            print(dt.today(), "bypassing:", message)
            return False
        else:
            print(dt.today(), "running:", message) 
            return True
    else:
        open(bypass_log_path, "a").close()
        print(dt.today(), "no bypass log.  running:", message)
        return True
        

            
        
# Used to determine quality encoding of fastq sequences.
# Assumes Phred+64 unless there is a character within the first 10000 reads with encoding in the Phred+33 range.
def check_code(segment):
    encoding = 64
    for item in segment:
        if(ord(item) < 64):
            encoding = 33
            break
    return encoding

def determine_encoding(fastq):
    #import the first 10k lines, then check the quality scores.
    #if the quality score symbols are below 76, it's phred33.  
    fastq_df = pd.read_csv(fastq, header=None, names=[None], sep="\n", skip_blank_lines = False, quoting=3, nrows=40000)
    fastq_df = pd.DataFrame(fastq_df.values.reshape(int(len(fastq_df)/4), 4))
    fastq_df.columns = ["ID", "seq", "junk", "quality"]
    quality_encoding = fastq_df["quality"].apply(lambda x: check_code(x)).mean() #condense into a single number.
    if(quality_encoding == 64): #all must be 64 or else it's 33
        quality_encoding = 64
    else:
        quality_encoding =  33
    return quality_encoding


# handles where to kill the pipeline, due to the prev step behaving badly
# logic is:  if the files inside the dep_path (or dep job label shortcut to the final_results)
#            are empty, then there's an error.  kill the pipeline 
def check_where_kill(dep_job_label=None, dep_path=None):
    if dep_job_label is None:
        if dep_path is None:
            return True
        else:
            dep_job_path = dep_path
    else:
        dep_job_path = os.path.join(dep_job_label, "final_results")

    file_list = os.listdir(dep_job_path)
    if len(file_list) > 0:
        for item in file_list:
            file_check_path = os.path.join(dep_job_path, item)
            if (os.path.getsize(file_check_path)) == 0:
                print("empty file detected: rerunning stage")
                sys.exit("bad dep")
        # run the job, silently
        return True
    else:
        print("stopping the pipeline.  dependencies don't exist")
        sys.exit("no dep")


# handles where to auto-resume the pipeline on a subsequent run
# label: used as a shorthand for paths we expect
# full path: a bypass for when we want to use it for detecting a location that doesn't fall into the normal format (final_results)
# dep: for checking if the job's dependencies are satisfied-> meant to point to the last stage's "final_results"
# logic is: if the full_path has no files (or the job label shortcut to final_results)
#           and the dependencies are ok, start the stage
#Aug 19, 2019: There's a tweak to this:  DIAMOND will generate zero-size files, due to no-matches
#it's allowable.
def check_where_resume(job_label=None, full_path=None, dep_job_path=None, file_check_bypass = False):
    if(not file_check_bypass):
        check_where_kill(dep_job_path)
    if job_label:
        job_path = os.path.join(job_label, "final_results")
    else:
        job_path = full_path

    print("looking at:", job_path)

    if os.path.exists(job_path):
        file_list = os.listdir(job_path)
        if(not file_check_bypass):
            if len(file_list) > 0:
                for item in file_list:
                    file_check_path = os.path.join(job_path, item)
                    if (os.path.getsize(file_check_path)) == 0:
                        print("empty file detected: rerunning stage")
                        return False
                print("bypassing!")
                return True
            else:
                print("no files detected: running")
                return False
        else:
            print("bypassing for special reasons")
            return True
    else:
        print("doesn't exist: running")
        return False




def main(config_path, pair_1_path, pair_2_path, single_path, output_folder_path, threads, args_pack):
    paths = mpp.tool_path_obj(config_path)
    no_host = args_pack["no_host"]
    keep_second_split = args_pack["keep_second_split"]
    verbose_mode = args_pack["verbose_mode"]
    infernal_limit = args_pack["infernal_limit"]
    chunks = paths.chunk_size
    BWA_mem_threshold = paths.BWA_mem_threshold
    BLAT_mem_threshold = paths.BLAT_mem_threshold
    DIAMOND_mem_threshold = paths.DIAMOND_mem_threshold
    
    
    if not single_path == "":
        read_mode = "single"
        quality_encoding = determine_encoding(single_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN SINGLE-ENDED MODE")
    else:
        read_mode = "paired"
        quality_encoding = determine_encoding(pair_1_path)
        print("ENCODING USED:", quality_encoding)
        print("OPERATING IN PAIRED-MODE")
        
    if threads == 0:
        real_thread_count = mp.cpu_count()
    else:
        real_thread_count = threads
       
    if(real_thread_count == 1):
        real_thread_count = 2
    print("number of threads used:", real_thread_count)         
            
    mp_store = []  # stores the multiprocessing processes

    # --------------------------------------------------
    # profiling vars
    # profiling vars are init here, in case a stage is skipped
    start_time = time.time()
    
    quality_start           = quality_end           = cleanup_quality_start             = cleanup_quality_end           = 0
    host_start              = host_end              = cleanup_host_start                = cleanup_host_end              = 0
    vector_start            = vector_end            = cleanup_vector_start              = cleanup_vector_end            = 0
    rRNA_filter_start       = rRNA_filter_end       = cleanup_rRNA_filter_start         = cleanup_rRNA_filter_end       = 0
    repop_start             = repop_end             = cleanup_repop_start               = cleanup_repop_end             = 0
    assemble_contigs_start  = assemble_contigs_end  = cleanup_assemble_contigs_start    = cleanup_assemble_contigs_end  = 0
    destroy_contigs_start   = destroy_contigs_end   = cleanup_destroy_contigs_start     = cleanup_destroy_contigs_end   = 0
    GA_BWA_start            = GA_BWA_end            = cleanup_GA_BWA_start              = cleanup_GA_BWA_end            = 0
    GA_BLAT_start           = GA_BLAT_end           = cleanup_GA_BLAT_start             = cleanup_GA_BLAT_end           = 0
    GA_DIAMOND_start        = GA_DIAMOND_end        = cleanup_GA_DIAMOND_start          = cleanup_GA_DIAMOND_end        = 0
    TA_start                = TA_end                = cleanup_TA_start                  = cleanup_TA_end                = 0
    EC_start                = EC_end                                                                                    = 0
    EC_DETECT_start         = EC_DETECT_end                                                                             = 0
    EC_PRIAM_start          = EC_PRIAM_end                                                                              = 0
    EC_DIAMOND_start        = EC_DIAMOND_end                                                                            = 0
    cleanup_EC_start        = cleanup_EC_end                                                                            = 0
    Cytoscape_start         = Cytoscape_end         = cleanup_cytoscape_start           = cleanup_cytoscape_end         = 0
    
    # the pipeline stages are all labelled.  This is for multiple reasons:  to keep the interim files organized properly
    # and to perform the auto-resume/kill features

    quality_filter_label                = "quality_filter"
    host_filter_label                   = "host_read_filter"
    vector_filter_label                 = "vector_read_filter"
    rRNA_filter_label                   = "rRNA_filter"
    rRNA_filter_second_split_label      = "rRNA_filter_second_split"
    rRNA_filter_barrnap_label           = "rRNA_filter_barrnap"
    rRNA_filter_infernal_label          = "rRNA_filter_infernal"
    rRNA_filter_post_label              = "rRNA_filter_post"
    repop_job_label                     = "duplicate_repopulation"
    assemble_contigs_label              = "assemble_contigs"
    destroy_contigs_label               = "destroy_contigs"
    gene_annotation_BWA_label           = "gene_annotation_BWA"
    gene_annotation_BWA_pp_label        = "gene_annotation_BWA_pp"
    gene_annotation_BLAT_label          = "gene_annotation_BLAT"
    gene_annotation_BLAT_cleanup_label  = "gene_annotation_BLAT_cleanup"
    gene_annotation_BLAT_cat_label      = "gene_annotation_BLAT_cat"
    gene_annotation_BLAT_pp_label       = "gene_annotation_BLAT_pp"
    gene_annotation_DIAMOND_label       = "gene_annotation_DIAMOND"
    gene_annotation_DIAMOND_pp_label    = "gene_annotation_DIAMOND_pp"
    gene_annotation_final_merge_label   = "gene_annotation_FINAL_MERGE"
    taxon_annotation_label              = "taxonomic_annotation"
    ec_annotation_label                 = "enzyme_annotation"
    ec_annotation_detect_label          = "enzyme_annotation_detect"
    ec_annotation_priam_label           = "enzyme_annotation_priam"
    ec_annotation_DIAMOND_label         = "enzyme_annotation_DIAMOND"
    ec_annotation_pp_label              = "enzyme_annotation_pp"
    output_label                        = "outputs"

    
    # Creates our command object, for creating shellscripts.
    if read_mode == "single":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, chunk_size = chunks, sequence_path_1=None, sequence_path_2=None, sequence_single=single_path)
    elif read_mode == "paired":
        commands = mpcom.mt_pipe_commands(no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, chunk_size = chunks, sequence_path_1=pair_1_path, sequence_path_2=pair_2_path, sequence_single=None)
    

    # This is the format we use to launch each stage of the pipeline.
    # We start a multiprocess that starts a subprocess.
    # The subprocess is created from the commands object

    # The quality filter stage
    quality_start = time.time()
    quality_path = os.path.join(output_folder_path, quality_filter_label)
    #if not check_where_resume(quality_path):
    if check_bypass_log(output_folder, quality_filter_label):
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                quality_filter_label,
                commands.create_quality_control_command(quality_filter_label),
                True
            )
        )
        process.start()  # start the multiprocess
        process.join()  # wait for it to end
        write_to_bypass_log(output_folder_path, quality_filter_label)
        cleanup_quality_start = time.time()
        if(verbose_mode == "quiet"):
            delete_folder(quality_path)
        elif(verbose_mode == "compress"):
            compress_folder(quality_path)
            delete_folder(quality_path)
        cleanup_quality_end = time.time()
    quality_end = time.time()
    print("quality filter:", '%1.1f' % (quality_end - quality_start - (cleanup_quality_end - cleanup_quality_start)), "s")
    print("quality filter cleanup:", '%1.1f' %(cleanup_quality_end - cleanup_quality_start), "s")
    

    # The host read filter stage
    if not no_host:
        host_start = time.time()
        host_path = os.path.join(output_folder_path, host_filter_label)
        #if not check_where_resume(host_path, None, quality_path):
        if check_bypass_log(output_folder_path, host_filter_label):
            process = mp.Process(
                target=commands.create_and_launch,
                args=(
                    host_filter_label,
                    commands.create_host_filter_command(host_filter_label, quality_filter_label),
                    True
                )
            )
            process.start()  # start the multiprocess
            process.join()  # wait for it to end
            write_to_bypass_log(output_folder_path, host_filter_label)
            cleanup_host_start = time.time()
            if(verbose_mode == "quiet"):
                delete_folder(host_path)
            elif(verbose_mode == "compress"):
                compress_folder(host_path)
                delete_folder(host_path)
            cleanup_host_end = time.time()
                
        host_end = time.time()
        print("host filter:", '%1.1f' % (host_end - host_start - (cleanup_host_end - cleanup_host_start)), "s")
        print("host filter cleanup:", '%1.1f' %(cleanup_host_end - cleanup_host_start),"s")
        

    # The vector contaminant filter stage
    vector_start = time.time()
    vector_path = os.path.join(output_folder_path, vector_filter_label)
    if no_host:
        #get dep args from quality filter
        #if not check_where_resume(vector_path, None, quality_path):
        if check_bypass_log(output_folder_path, host_filter_label):
            process = mp.Process(
                target=commands.create_and_launch,
                args=(
                    vector_filter_label,
                    commands.create_vector_filter_command(vector_filter_label, quality_filter_label),
                    True
                )
            )
            process.start()  # start the multiprocess
            process.join()  # wait for it to end
            write_to_bypass_log(output_folder_path, vector_filter_label)
            cleanup_vector_start = time.time()
            if(verbose_mode == "quiet"):
                delete_folder(vector_path)
            elif(verbose_mode == "compress"):
                compress_folder(vector_path)
                delete_folder(vector_path)
            cleanup_vector_end = time.time()
    else:
        #get the dep args from host filter
        #if not check_where_resume(vector_path, None, host_path):
        if check_bypass_log(output_folder_path, vector_filter_label):
            process = mp.Process(
                target=commands.create_and_launch,
                args=(
                    vector_filter_label,
                    commands.create_vector_filter_command(vector_filter_label, host_filter_label),
                    True
                )
            )
            process.start()  # start the multiprocess
            process.join()  # wait for it to end
            write_to_bypass_log(output_folder_path, vector_filter_label)
            cleanup_vector_start = time.time()
            if(verbose_mode == "quiet"):
                delete_folder(vector_path)
            elif(verbose_mode == "compress"):
                compress_folder(vector_path)
                delete_folder(vector_path)
            cleanup_vector_end = time.time()
    vector_end = time.time()
    print("vector filter:", '%1.1f' % (vector_end - vector_start - (cleanup_vector_end - cleanup_vector_start)), "s")
    print("vector filter cleanup:", '%1.1f' % (cleanup_vector_end - cleanup_vector_start), "s")
    

    # ----------------------------------------------
    # rRNA removal stage
    rRNA_filter_start = time.time()

    rRNA_filter_path = os.path.join(output_folder_path, rRNA_filter_label)
    #if not check_where_resume(rRNA_filter_path, None, vector_path):
    if check_bypass_log(output_folder_path, rRNA_filter_label): 
        sections = ["singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in reversed(sections):  #we go backwards due to a request by Ana.  pairs first, if applicable, then singletons
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            split_path = os.path.join(rRNA_filter_path, "data", section, section + "_fastq")
            second_split_path = os.path.join(rRNA_filter_path, "data", section, section + "_second_split_fastq")
            barrnap_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section, section + "_barrnap")
            infernal_path = os.path.join(output_folder_path, rRNA_filter_label, "data", section, section + "_infernal") 
            
            #if not check_where_resume(job_label = None, full_path = second_split_path, dep_job_path = vector_path):
            if check_bypass_log(output_folder_path, rRNA_filter_second_split_label + "_" + section):
                print(dt.today(), "splitting:", section, " for rRNA filtration")
                inner_name = "rRNA_filter_prep_" + section
                process = mp.Process(
                    target = commands.create_and_launch,
                    args = (
                        rRNA_filter_label,
                        commands.create_rRNA_filter_prep_command_v2(rRNA_filter_label, section, vector_filter_label),
                        True,
                        inner_name
                    )
                )
        
                process.start()
                process.join()
                write_to_bypass_log(output_folder_path, rRNA_filter_second_split_label + "_" + section)
                
            #secondary split -> number of files
                concurrent_job_count = 0
                cumulative_job_count = 0
                batch_count = 0
                
                for item in os.listdir(split_path):
                    #second_split_path = os.path.split(
                    inner_name = "rRNA_filter_prep_2_" + section + "_" + str(cumulative_job_count)
                    full_item_path = os.path.join(split_path, item)
                    print(dt.today(), "full item path:", full_item_path)
                    cumulative_job_count += 1
                    process = mp.Process(
                        target = commands.create_and_launch,
                        args = (
                            rRNA_filter_label,
                            commands.create_rRNA_filter_prep_command_2nd_split(rRNA_filter_label, section, full_item_path,  40),
                            True,
                            inner_name
                        )
                    )
                    mp_store.append(process)
                    process.start()
                    concurrent_job_count += 1
                    if(concurrent_job_count >= infernal_limit):
                        for p_item in mp_store:
                            p_item.join()
                        mp_store[:] = []
                        concurrent_job_count = 0
                        print(dt.today(), "waiting for rRNA second split to finish running.  batch:", batch_count)
                        batch_count += 1
                        print(dt.today(), "pausing for 2 seconds to let things settle")
                        time.sleep(2)
                        print(dt.today(), "pause over.  resuming")
                        
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []
        
        
            #if not check_where_resume(job_label = None, full_path = barrnap_path, dep_job_path = vector_path):
            if check_bypass_log(output_folder_path, rRNA_filter_barrnap_label + "_" + section):
                concurrent_job_count = 0
                batch_count = 0
                for item in os.listdir(second_split_path):
                    print(dt.today(), "barrnap looking at:", item)
                    inner_name = "rRNA_filter_barrnap_" + item.split(".")[0]
                    print(dt.today(), "barrnap job script name", inner_name)
                    print("rRNA filter inner name:", inner_name)
                    
                    process = mp.Process(
                        target=commands.create_and_launch,
                        args=(
                            rRNA_filter_label,
                            commands.create_rRNA_filter_barrnap_command("rRNA_filter", section, item, vector_filter_label),
                            True,
                            inner_name
                        )
                    )
                    mp_store.append(process)
                    process.start()
                    
                    concurrent_job_count += 1
                    if(concurrent_job_count >= infernal_limit): 
                        print(dt.today(), "letting a batch job run: barrnap", batch_count)
                        for p_item in mp_store:
                            p_item.join()
                        mp_store[:] = []  # clear the list    
                        concurrent_job_count = 0
                        batch_count += 1
                        
                        print(dt.today(), "pausing barrnap for 2 seconds to let things settle")
                        time.sleep(2)
                        print(dt.today(), "brakes off.  continue barrnap")
                        
                print(dt.today(), "final batch: barrnap")
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []  # clear the list    
                write_to_bypass_log(output_folder_path, rRNA_filter_barrnap_label + "_" + section)
                
                
            #if not check_where_resume(job_label = None, full_path = infernal_path):#, dep_job_path = barrnap_path):
            if check_bypass_log(output_folder_path, rRNA_filter_infernal_label + "_" + section):
                concurrent_job_count = 0
                batch_count = 0
                #these jobs now have to be launched in segments
                for item in os.listdir(second_split_path):
                    inner_name = "rRNA_filter_infernal_" + item.split(".")[0]
                    print("rRNA filter inner name:", inner_name)
                    
                    process = mp.Process(
                        target=commands.create_and_launch,
                        args=(
                            rRNA_filter_label,
                            commands.create_rRNA_filter_infernal_command("rRNA_filter", section, item, vector_filter_label),
                            True,
                            inner_name
                        )
                    )
                    mp_store.append(process)
                    process.start()
                    concurrent_job_count += 1
                    if(concurrent_job_count >= infernal_limit):
                        
                        print(dt.today(), "letting small batch of jobs run: infernal", batch_count)
                        for p_item in mp_store:
                            p_item.join()
                        mp_store[:] = []  # clear the list
                        batch_count += 1
                        concurrent_job_count = 0
                        
                print(dt.today(), "final batch: infernal")
                for p_item in mp_store:
                    p_item.join()
                mp_store[:] = []  # clear the list
                write_to_bypass_log(output_folder_path, rRNA_filter_infernal_label + "_" + section)
            
        if check_bypass_log(output_folder_path, rRNA_filter_post_label):
            print(dt.today(), "now running rRNA filter post")
            inner_name = "rRNA_filter_post"
            process = mp.Process(
                target=commands.create_and_launch,
                args=(
                    rRNA_filter_label,
                    commands.create_rRNA_filter_post_command(vector_filter_label, rRNA_filter_label),
                    True,
                    inner_name
                )
            )
            process.start()
            process.join()
            write_to_bypass_log(output_folder_path, rRNA_filter_post_label)
        
        if not (keep_second_split):
            print("deleting rRNA second split")
            for item in sections:
                second_split_path = os.path.join(rRNA_filter_path, "data", section, section + "_second_split_fastq")
                if(os.path.exists(second_split_path)):
                    shutil.rmtree(second_split_path)
                else:
                    print(dt.today(), "already deleted:", second_split_path)
        write_to_bypass_log(output_folder_path, rRNA_filter_label)
        cleanup_rRNA_filter_start = time.time()
        if(verbose_mode == "quiet"):
            delete_folder(rRNA_filter_path)
        elif(verbose_mode == "compress"):
            compress_folder(rRNA_filter_path)
            delete_folder(rRNA_filter_path)
        cleanup_rRNA_filter_end = time.time()
    rRNA_filter_end = time.time()
    
    print("rRNA filter:", '%1.1f' % (rRNA_filter_end - rRNA_filter_start - (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start)), "s")
    print("rRNA filter cleanup:", '%1.1f' % (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start), "s")
    
    # Duplicate repopulation
    repop_start = time.time()
    repop_job_path = os.path.join(output_folder_path, repop_job_label)
    #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
    if check_bypass_log(output_folder_path, repop_job_label):
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                repop_job_label,
                commands.create_repop_command(repop_job_label, quality_filter_label, rRNA_filter_label),
                True
            )
        )
        process.start()
        process.join()
        write_to_bypass_log(output_folder_path, repop_job_label)
        
    cleanup_repop_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(repop_job_path)
    elif(verbose_mode == "compress"):
        compress_folder(repop_job_path)
        delete_folder(repop_job_path)
    cleanup_repop_end = time.time()
    repop_end = time.time()
    print("repop:", '%1.1f' % (repop_end - repop_start - (cleanup_repop_end - cleanup_repop_start)), "s")
    print("repop cleanup:", '%1.1f' % (cleanup_repop_end - cleanup_repop_start), "s")
    # -------------------------------------------------------------
    
    # ----------------------------------------
    # Assemble contigs
    assemble_contigs_start = time.time()
    assemble_contigs_path = os.path.join(output_folder_path, assemble_contigs_label)
    
    
    #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
    if check_bypass_log(output_folder_path, assemble_contigs_label):
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                assemble_contigs_label,
                commands.create_assemble_contigs_command(assemble_contigs_label, repop_job_label),
                True
            )
        )
        process.start()
        process.join()
        
        mgm_file = os.path.join(assemble_contigs_path, "data", "1_mgm", "gene_report.txt")
        if(os.path.exists(mgm_file)):
            write_to_bypass_log(output_folder_path, assemble_contigs_label)
        else:
            sys.exit("mgm did not run.  look into it.  pipeline stopping here")
        
        
        
       
    cleanup_assemble_contigs_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(assemble_contigs_path)
    elif(verbose_mode == "compress"):
        compress_folder(assemble_contigs_path)
        delete_folder(assemble_contigs_path)
    cleanup_assemble_contigs_end = time.time()
    assemble_contigs_end = time.time()
    print("assemble contigs:", '%1.1f' % (assemble_contigs_end - assemble_contigs_start - (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start)), "s")    
    print("assemble contigs cleanup:", '%1.1f' % (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start), "s")

    
    # ----------------------------------------------
    # BWA gene annotation
    
    
    
    
    GA_BWA_start = time.time()
    gene_annotation_BWA_path = os.path.join(output_folder_path, gene_annotation_BWA_label)
    #if not check_where_resume(gene_annotation_BWA_path, None, assemble_contigs_path):
    if check_bypass_log(output_folder_path, gene_annotation_BWA_label):
        process = mp.Process(
            target = commands.create_and_launch, 
            args = (
                gene_annotation_BWA_label, 
                commands.create_split_ga_fasta_data_command(gene_annotation_BWA_label, assemble_contigs_label, "contigs"),
                True, 
                "GA_prep_split_contigs"
            )
        )
        process.start()
        mp_store.append(process)
        sections = ["singletons"]
        if(read_mode == "paired"):
            sections.extend(["pair_1", "pair_2"])
        for section in sections:    
            process = mp.Process(
                target = commands.create_and_launch, 
                args = (
                    gene_annotation_BWA_label, 
                    commands.create_split_ga_fastq_data_command(gene_annotation_BWA_label, assemble_contigs_label, section),
                    True,
                    "GA_prep_split_" + section
                )
            )
            process.start()
            mp_store.append(process)
        
        for item in mp_store:
            item.join()
        mp_store[:] = []
        
        #-------------------------------------------------------------------------
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(gene_annotation_BWA_path, "data", "0_read_split", section)):
                job_submitted = False
                full_sample_path = os.path.join(os.path.join(gene_annotation_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA" + "_" + file_tag
                BWA_output_name = file_tag + ".sam"
                BWA_output_path = os.path.join(gene_annotation_BWA_path, "data", "1_bwa", BWA_output_name)
                #this checker assumes that BWA only exports a file when it's finished running
                if(os.path.exists(BWA_output_path)):
                    continue
                else:
                    while(not job_submitted): 
                        if(mem_checker(BWA_mem_threshold)):
                            print(dt.today(), "mem ok:", psu.virtual_memory().available / (1024*1024*1000), "GB")
                            
                            bwa_process = mp.Process(
                                target = commands.create_and_launch, 
                                args = (gene_annotation_BWA_label, commands.create_BWA_annotate_command_v2(gene_annotation_BWA_label, full_sample_path), True, job_name)
                            )
                            bwa_process.start()
                            
                            mp_store.append(bwa_process)
                            job_submitted = True
                            print(dt.today(), "submitted:", job_name)
                            time.sleep(5) #placed here so the process has some time to get started.
                            
                        else:               
                            print(dt.today(), "mem has reached a limit.  waiting:", job_name, "available mem:", psu.virtual_memory().available / (1024*1024*1000), "GB")
                            
                            time.sleep(5)
                            #for item in mp_store:
                            #    item.join()
                            #mp_store[:] = []
        print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
        for item in mp_store:
            item.join()
        mp_store[:] = []

       
        
        write_to_bypass_log(output_folder_path, gene_annotation_BWA_label)
        
    
    if check_bypass_log(output_folder_path, gene_annotation_BWA_pp_label):
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(gene_annotation_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(gene_annotation_BWA_path, "data", "0_read_split", section, split_sample))
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA_pp" + "_" + file_tag
                job_submitted = False
                while(not job_submitted):
                    if(mem_checker(10)):
                        print(dt.today(), "BWA PP mem ok:", psu.virtual_memory().available / (1024*1024*1000), "GB")
                        BWA_pp_process = mp.Process(target = commands.create_and_launch,
                            args = (
                                gene_annotation_BWA_label, 
                                commands.create_BWA_pp_command_v2(gene_annotation_BWA_label, assemble_contigs_label, full_sample_path),
                                True,
                                job_name
                            )
                        )
                        BWA_pp_process.start()
                        mp_store.append(BWA_pp_process)
                        job_submitted = True
                        print(dt.today(), "submitted:", job_name)
                        time.sleep(1)
                    else:
                        time.sleep(5)
                        print(dt.today(), "BWA pp mem busy:", psu.virtual_memory().available / (1024*1024*1000), "GB")
                        
                        
        print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
        for item in mp_store:
            item.join()
        mp_store[:] = []

        
        process = mp.Process(
            target = commands.create_and_launch,
            args = (
                gene_annotation_BWA_label, commands.create_BWA_copy_contig_map_command(gene_annotation_BWA_label, assemble_contigs_label),
                True,
                "GA_copy_contigs"
            )
        )
        process.start()
        process.join()
        
        
        write_to_bypass_log(output_folder_path, gene_annotation_BWA_pp_label)
        
    cleanup_GA_BWA_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(gene_annotation_BWA_path)
    elif(verbose_mode == "compress"):
        compress_folder(gene_annotation_BWA_path)
        delete_folder(gene_annotation_BWA_path)
    cleanup_GA_BWA_end = time.time()
    GA_BWA_end = time.time()
    print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
    print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
    
    # ------------------------------------------------
    # BLAT gene annotation
    GA_BLAT_start = time.time()
    gene_annotation_BLAT_path = os.path.join(output_folder_path, gene_annotation_BLAT_label)
    #if not check_where_resume(gene_annotation_BLAT_path, None, gene_annotation_BWA_path):
    if check_bypass_log(output_folder_path, gene_annotation_BLAT_label):
        BlatPool = mp.Pool(int(real_thread_count / 2))
        print(dt.today(), "BLAT threads used:", real_thread_count/2)
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        sample_job_flag = True
        missed_jobs_list = []
        
        for section in sections:
            for split_sample in os.listdir(os.path.join(gene_annotation_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(gene_annotation_BWA_path, "final_results", split_sample))
                    for fasta_db in os.listdir(paths.DNA_DB_Split):
                        if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                            job_name = "BLAT_" + file_tag + "_" + fasta_db
                            blat_output_name = file_tag + "_" + fasta_db + ".blatout"
                            blat_output_path = os.path.join(gene_annotation_BLAT_path, "data", "0_blat", blat_output_name)
                            #This checker assume BLAT only exports a file when it's finished running
                            if(os.path.exists(blat_output_path)):
                                #print(dt.today(), "BLAT job ran already, skipping:", blat_output_name)
                                continue
                            else:
                                job_submitted = False
                                while(not job_submitted):
                                    if(mem_checker(BLAT_mem_threshold)):
                                        #why this? because we have too many BLAT jobs, and the deletion of those shellscript files added too much of a delay.
                                        BLAT_process = mp.Process(target = commands.launch_only,#commands.create_and_launch,
                                            args=(
                                                gene_annotation_BLAT_label,
                                                commands.create_BLAT_annotate_command_v2(gene_annotation_BLAT_label, full_sample_path, fasta_db),
                                                True,
                                                job_name
                                            )
                                        )
                                        BLAT_process.start()
                                        mp_store.append(BLAT_process)
                                        job_submitted = True
                                        print(dt.today(), job_name, "submitted: mem ok:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                                    else:
                                        print(dt.today(), "too many BLAT jobs. waiting:", job_name, psu.virtual_memory().available/(1024*1024*1000), "GB")
                                        time.sleep(5)
                                
        print(dt.today(), "final BLAT job removal")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        #for item in os.listdir(gene_annotation_BLAT_path):
        #    if(item.endswith(".ffn.sh")):
        #        if(os.path.exists(item)):
        #            os.remove(item)

        
        write_to_bypass_log(output_folder_path, gene_annotation_BLAT_label)
        
    #retired    
    if check_bypass_log(output_folder_path, gene_annotation_BLAT_cleanup_label):
        for item in os.listdir(gene_annotation_BLAT_path):
            if(item.endswith(".ffn.sh")):
                job_path = os.path.join(gene_annotation_BLAT_path, item)
                os.remove(job_path)
        write_to_bypass_log(output_folder_path, gene_annotation_BLAT_cleanup_label)
        
    if check_bypass_log(output_folder_path, gene_annotation_BLAT_cat_label):
        for split_sample in os.listdir(os.path.join(gene_annotation_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                full_sample_path = os.path.join(os.path.join(gene_annotation_BWA_path, "final_results", split_sample))
                job_name = file_tag + "_cat"
                process = mp.Process(
                    target=commands.create_and_launch,
                    args=(
                        gene_annotation_BLAT_label,
                        commands.create_BLAT_cat_command_v2(gene_annotation_BLAT_label, full_sample_path),
                        True,
                        job_name
                    )
                )
                process.start()
                mp_store.append(process)
                
            
        for item in mp_store:
            item.join()
        mp_store[:] = []
        write_to_bypass_log(output_folder_path, gene_annotation_BLAT_cat_label)
    
    if check_bypass_log(output_folder_path, gene_annotation_BLAT_pp_label):
        for split_sample in os.listdir(os.path.join(gene_annotation_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BLAT_" + file_tag + "_pp"
                full_sample_path = os.path.join(os.path.join(gene_annotation_BWA_path, "final_results", split_sample))
                job_submitted = False
                while(not job_submitted):
                    if(mem_checker(10)):
                        Blat_pp_process = mp.Process(target = commands.create_and_launch,
                            args=(
                                gene_annotation_BLAT_label,
                                commands.create_BLAT_pp_command_v2(gene_annotation_BLAT_label, full_sample_path, gene_annotation_BWA_label),
                                True,
                                job_name
                            )
                        )
                        Blat_pp_process.start()
                        mp_store.append(Blat_pp_process)
                        print(dt.today(), job_name, "submitted. mem:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                        time.sleep(1)
                        job_submitted = True
                    else:
                        print(dt.today(), job_name, "waiting. mem:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                        time.sleep(5)
                        
        print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
 
        process = mp.Process(
            target = commands.create_and_launch,
            args = (
                gene_annotation_BLAT_label, commands.create_BLAT_copy_contig_map_command(gene_annotation_BLAT_label, gene_annotation_BWA_label),
                True,
                "GA_BLAT_copy_contigs"
            )
        )
        process.start()
        process.join()
        
        write_to_bypass_log(output_folder_path, gene_annotation_BLAT_pp_label)

    cleanup_GA_BLAT_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(gene_annotation_BLAT_path)
    elif(verbose_mode == "compress"):
        compress_folder(gene_annotation_BLAT_path)
        delete_folder(gene_annotation_BLAT_path)
    cleanup_GA_BLAT_end = time.time()
    GA_BLAT_end = time.time()
    print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
    print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
    
    # ------------------------------------------------------
    # Diamond gene annotation
    GA_DIAMOND_start = time.time()
    gene_annotation_DIAMOND_path = os.path.join(output_folder_path, gene_annotation_DIAMOND_label)
    GA_DIAMOND_tool_output_path = os.path.join(gene_annotation_DIAMOND_path, "data", "0_diamond")
    #if not check_where_resume(None, GA_DIAMOND_tool_output_path, gene_annotation_BLAT_path, file_check_bypass = True):
    if check_bypass_log(output_folder_path, gene_annotation_DIAMOND_label):
        for split_sample in os.listdir(os.path.join(gene_annotation_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_" + file_tag
                full_sample_path = os.path.join(os.path.join(gene_annotation_BLAT_path, "final_results", split_sample))
                job_submitted = False
                while(not job_submitted):
                    if(mem_checker(DIAMOND_mem_threshold)):
                        DIAMOND_process = mp.Process(target = commands.create_and_launch,
                            args=(
                                gene_annotation_DIAMOND_label,
                                commands.create_DIAMOND_annotate_command_v2(gene_annotation_DIAMOND_label, full_sample_path),
                                True,
                                job_name
                            )
                        )
                        DIAMOND_process.start()
                        mp_store.append(DIAMOND_process)
                        print(dt.today(), "Submitted:", job_name, psu.virtual_memory().available/(1024*1024*1000), "GB")
                        job_submitted = True
                        time.sleep(5) #it's enough time for the process to eat some memory.
                        
                    else:
                        print(dt.today(), "DIAMOND: we've reached the mem limit. waiting:", job_name)
                        time.sleep(5)
                        #for item in mp_store:
                        #    item.join()
                        #mp_store[:] = []
        
        # DIAMOND_Pool = mp.Pool(int(real_thread_count / 2))
        # print(dt.today(), "DIAMOND threads used:", real_thread_count/2)
        # for split_sample in os.listdir(os.path.join(gene_annotation_BLAT_path, "final_results")):
            # if(split_sample.endswith(".fasta")):
                # file_tag = os.path.basename(split_sample)
                # file_tag = os.path.splitext(file_tag)[0]
                # job_name = "DIAMOND_" + file_tag
                # full_sample_path = os.path.join(os.path.join(gene_annotation_BLAT_path, "final_results", split_sample))
                # DIAMOND_Pool.apply_async(commands.create_and_launch,
                    # args=(
                        # gene_annotation_DIAMOND_label,
                        # commands.create_DIAMOND_annotate_command_v2(gene_annotation_DIAMOND_label, full_sample_path),
                        # True,
                        # job_name
                    # )
                # )
        # DIAMOND_Pool.close()
        # DIAMOND_Pool.join()
        print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
        for item in mp_store:
            item.join()
        mp_store[:] = []
        write_to_bypass_log(output_folder_path, gene_annotation_DIAMOND_label)
        
    #if not check_where_resume(gene_annotation_DIAMOND_path, None, GA_DIAMOND_tool_output_path, file_check_bypass = True):
    if check_bypass_log(output_folder_path, gene_annotation_DIAMOND_pp_label):
        DIAMOND_pp_Pool = mp.Pool(int(real_thread_count / 2))
        print(dt.today(), "DIAMOND PP threads used:", real_thread_count/2)
        for split_sample in os.listdir(os.path.join(gene_annotation_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_pp_" + file_tag
                full_sample_path = os.path.join(os.path.join(gene_annotation_BLAT_path, "final_results", split_sample))
                job_submitted = False
                while(not job_submitted):
                    if(mem_checker(10)):
                        DIAMOND_pp_process = mp.Process(target = commands.create_and_launch,
                            args=(
                                gene_annotation_DIAMOND_label,
                                commands.create_DIAMOND_pp_command_v2(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label, full_sample_path),
                                True,
                                job_name
                            )
                        )
                        DIAMOND_pp_process.start()
                        mp_store.append(DIAMOND_pp_process)
                        job_submitted = True
                        print(dt.today(), "Submitted:", job_name, psu.virtual_memory().available/(1024*1024*1000), "GB")
                        time.sleep(1)
                    else:
                        
                        print(dt.today(), job_name, "waiting.  mem:", psu.virtual_memory().available/(1024*1024*1000), "GB")
                        time.sleep(5)
                        
        print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
        for item in mp_store:
            item.join()
        mp_store[:] = []
                # DIAMOND_pp_Pool.apply_async(commands.create_and_launch,
                    # args=(
                        # gene_annotation_DIAMOND_label,
                        # commands.create_DIAMOND_pp_command_v2(gene_annotation_DIAMOND_label, gene_annotation_BLAT_label, full_sample_path),
                        # True,
                        # job_name
                    # )
                # )
        # DIAMOND_pp_Pool.close()
        # DIAMOND_pp_Pool.join()
            
        write_to_bypass_log(output_folder_path, gene_annotation_DIAMOND_pp_label)
    
        
    
    cleanup_GA_DIAMOND_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(gene_annotation_DIAMOND_path)
    elif(verbose_mode == "compress"):
        compress_folder(gene_annotation_DIAMOND_path)
        delete_folder(gene_annotation_DIAMOND_path)
    cleanup_GA_DIAMOND_end = time.time()
    GA_DIAMOND_end = time.time()
    print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
    print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
    
    
    GA_final_merge_start = time.time()
    if check_bypass_log(output_folder_path, gene_annotation_final_merge_label):
        final_merge_process = mp.Process(
            target = commands.create_and_launch, 
            args = (
                gene_annotation_final_merge_label,
                commands.create_GA_final_merge_command(gene_annotation_final_merge_label, gene_annotation_BWA_label, gene_annotation_BLAT_label, gene_annotation_DIAMOND_label, assemble_contigs_label),
                True, 
                "GA_final_merge"
            )
        )
        final_merge_process.start()
        final_merge_process.join()
        write_to_bypass_log(output_folder_path, gene_annotation_final_merge_label)
    GA_final_merge_end = time.time()
    print("GA final merge:", '%1.1f' % (GA_final_merge_end - GA_final_merge_start), "s")
    
    # ------------------------------------------------------
    # Taxonomic annotation
    TA_start = time.time()
    taxon_annotation_path = os.path.join(output_folder_path, taxon_annotation_label)
    #if not check_where_resume(taxon_annotation_path, None, gene_annotation_DIAMOND_path):
    if check_bypass_log(output_folder_path, taxon_annotation_label):
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                taxon_annotation_label,
                commands.create_taxonomic_annotation_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, gene_annotation_final_merge_label),
                True
            )
        )
        process.start()
        process.join()
        write_to_bypass_log(output_folder_path, taxon_annotation_label)
    cleanup_TA_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(taxon_annotation_path)
    elif(verbose_mode == "compress"):
        compress_folder(taxon_annotation_path)
        delete_folder(taxon_annotation_path)
    cleanup_TA_end = time.time()
    TA_end = time.time()
    print("TA:", '%1.1f' % (TA_end - TA_start - (cleanup_TA_end - cleanup_TA_start)), "s")
    print("TA cleanup:", '%1.1f' % (cleanup_TA_end - cleanup_TA_start), "s")
    
    
    
    
    
    
    
    # ------------------------------------------------------
    # Detect EC annotation
    EC_process_list = []
    
    ec_annotation_path = os.path.join(output_folder_path, ec_annotation_label)
    EC_start = time.time()
    #There's a 2-step check.  We don't want it ti re-run either DETECT, or PRIAM+DIAMOND because they're too slow
    #if not check_where_resume(ec_annotation_path, None, gene_annotation_DIAMOND_path):
    #if check_bypass_log(output_folder_path, ec_annotation_label):
    EC_DETECT_start = time.time()
    ec_detect_path = os.path.join(ec_annotation_path, "data", "0_detect")
    #if not check_where_resume(job_label = None, full_path = ec_detect_path, dep_job_path = gene_annotation_DIAMOND_path):
    if check_bypass_log(output_folder_path, ec_annotation_detect_label):
        inner_name = "ec_detect"
        process = mp.Process(
            target = commands.create_and_launch,
            args = (
                ec_annotation_label, 
                commands.create_EC_DETECT_command(ec_annotation_label, gene_annotation_final_merge_label),
                True,
                inner_name
            )
        )
        process.start()
        EC_process_list.append(process)
        #process.join()
        
    EC_DETECT_end = time.time()
    print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    
    # --------------------------------------------------------------
    # Priam EC annotation.  Why isn't it parallel? computing restraints.  Not enough mem
    EC_PRIAM_start = time.time()
    
    ec_priam_path = os.path.join(ec_annotation_path, "data", "1_priam")
    #if not check_where_resume(job_label = None, full_path = ec_priam_path, dep_job_path = gene_annotation_DIAMOND_path):
    if check_bypass_log(output_folder_path, ec_annotation_priam_label):
        inner_name = "ec_priam"
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                ec_annotation_label,
                commands.create_EC_PRIAM_command(ec_annotation_label, gene_annotation_final_merge_label),
                True,
                inner_name
            )
        )
        process.start()
        EC_process_list.append(process)
        #process.join()
    EC_PRIAM_end = time.time()
    print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
    # --------------------------------------------------------------
    # DIAMOND EC annotation 
    EC_DIAMOND_start = time.time()
    ec_diamond_path = os.path.join(ec_annotation_path, "data", "2_diamond")
    #if not check_where_resume(job_label = None, full_path = ec_diamond_path, dep_job_path = gene_annotation_DIAMOND_path):
    if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        inner_name = "ec_diamond"
        process = mp.Process(
            target = commands.create_and_launch, 
            args = (
                ec_annotation_label,
                commands.create_EC_DIAMOND_command(ec_annotation_label, gene_annotation_final_merge_label),
                True,
                inner_name
            )
        )
        process.start()
        EC_process_list.append(process)
        #process.join()
    EC_DIAMOND_end = time.time()
    for item in EC_process_list:
        item.join()
    EC_process_list[:] = []
    
    ec_detect_out = os.path.join(ec_detect_path, "proteins.fbeta")
    ec_priam_out = os.path.join(ec_priam_path, "PRIAM_proteins_priam", "ANNOTATION", "sequenceECs.txt")
    ec_diamond_out = os.path.join(ec_diamond_path, "proteins.blastout")
    if check_bypass_log(output_folder_path, ec_annotation_detect_label):
        if(os.path.exists(ec_detect_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_detect_label)
    if check_bypass_log(output_folder_path, ec_annotation_priam_label):
        if(os.path.exists(ec_priam_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_priam_label)
    if check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        if(os.path.exists(ec_diamond_out)):
            write_to_bypass_log(output_folder_path, ec_annotation_DIAMOND_label)
    
    #----------------------------------------------------------------------
    # EC post process
    EC_post_start = time.time()
    #if not (check_where_resume(ec_annotation_path, None, gene_annotation_DIAMOND_path)):
    if check_bypass_log(output_folder_path, ec_annotation_pp_label):
        inner_name = "ec_post"
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                ec_annotation_label,
                commands.create_EC_postprocess_command(ec_annotation_label, gene_annotation_final_merge_label),
                True,
                inner_name
            )
        )
        process.start()
        process.join()
        write_to_bypass_log(output_folder_path, ec_annotation_pp_label)
    cleanup_EC_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(ec_annotation_path)
    elif(verbose_mode == "compress"):
        compress_folder(ec_annotation_path)
        delete_folder(ec_annotation_path)
    cleanup_EC_end = time.time()
    EC_post_end = time.time()
        
    #else:
    #    #EC bypassed
    #    EC_PRIAM_start = time.time()
    #    EC_PRIAM_end = time.time()
    #    EC_DIAMOND_start = time.time()
    #    EC_DIAMOND_end = time.time()
    #    EC_DETECT_start = time.time()
    #    EC_DETECT_end = time.time()
    #    EC_post_start = time.time()
    #    EC_post_end = time.time()
    #    cleanup_EC_start = time.time()
    #    cleanup_EC_end = time.time()
        
        #print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    #EC_PRIAM_DIAMOND_end = time.time()
    #print("EC PRIAM + DIAMOND:", '%1.1f' % (EC_PRIAM_DIAMOND_end - EC_PRIAM_DIAMOND_start - (cleanup_EC_end - cleanup_EC_start)), "s")
    EC_end = time.time()
    print("EC run:", '%1.1f' % (EC_end - EC_start), "s")
    print("EC cleanup:", '%1.1f' % (cleanup_EC_end - cleanup_EC_start), "s")
    
    # ------------------------------------------------------
    # RPKM Table and Cytoscape Network
    Cytoscape_start = time.time()
    network_path = os.path.join(output_folder_path, output_label)
    #if not check_where_resume(network_path, None, ec_annotation_path):
    if check_bypass_log(output_folder_path, output_label):
        process = mp.Process(
            target=commands.create_and_launch,
            args=(
                output_label,
                commands.create_output_generation_command(output_label, quality_filter_label, host_filter_label, 
                assemble_contigs_label, repop_job_label, gene_annotation_final_merge_label, taxon_annotation_label, ec_annotation_label), 
                True
            )
        )
        process.start()
        process.join()
        write_to_bypass_log(output_folder_path, output_label)
    cleanup_cytoscape_start = time.time()
    if(verbose_mode == "quiet"):
        delete_folder(network_path)
    elif(verbose_mode == "compress"):
        compress_folder(network_path)
        delete_folder(network_path)
    cleanup_cytoscape_end = time.time()
        
    Cytoscape_end = time.time()
    end_time = time.time()
    print("Outputs:", '%1.1f' % (Cytoscape_end - Cytoscape_start - (cleanup_cytoscape_end - cleanup_cytoscape_start)), "s")
    print("Outputs cleanup:", '%1.1f' % (cleanup_cytoscape_end - cleanup_cytoscape_start), "s")
    print("=============================================================================================")
    print("Final summary")
    print("--------------------------------------------------------")
    print("Total runtime:", '%1.1f' % (end_time - start_time), "s")
    print("quality filter:", '%1.1f' % (quality_end - quality_start - (cleanup_quality_end - cleanup_quality_start)), "s")
    print("quality filter cleanup:", '%1.1f' %(cleanup_quality_end - cleanup_quality_start), "s")
    if not no_host:
        print("host filter:", '%1.1f' % (host_end - host_start - (cleanup_host_end - cleanup_host_start)), "s")
        print("host filter cleanup:", '%1.1f' %(cleanup_host_end - cleanup_host_start),"s")
    print("vector filter:", '%1.1f' % (vector_end - vector_start - (cleanup_vector_end - cleanup_vector_start)), "s")
    print("vector filter cleanup:", '%1.1f' % (cleanup_vector_end - cleanup_vector_start), "s")
    print("rRNA filter:", '%1.1f' % (rRNA_filter_end - rRNA_filter_start - (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start)), "s")
    print("rRNA filter cleanup:", '%1.1f' % (cleanup_rRNA_filter_end - cleanup_rRNA_filter_start), "s")
    print("repop:", '%1.1f' % (repop_end - repop_start - (cleanup_repop_end - cleanup_repop_start)), "s")
    print("repop cleanup:", '%1.1f' % (cleanup_repop_end - cleanup_repop_start), "s")
    print("assemble contigs:", '%1.1f' % (assemble_contigs_end - assemble_contigs_start - (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start)), "s")    
    print("assemble contigs cleanup:", '%1.1f' % (cleanup_assemble_contigs_end - cleanup_assemble_contigs_start), "s")
    print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
    print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
    print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
    print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
    print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
    print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
    print("TA:", '%1.1f' % (TA_end - TA_start - (cleanup_TA_end - cleanup_TA_start)), "s")
    print("TA cleanup:", '%1.1f' % (cleanup_TA_end - cleanup_TA_start), "s")
    print("EC:", '%1.1f' % (EC_end - EC_start), "s")
    #print("---------------------------------------------")
    #print("Note: EC is in cloud-mode.  ignore individual timing")
    #print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    #print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
    #print("EC DIAMOND:", '%1.1f' % (EC_DIAMOND_end - EC_DIAMOND_start), "s")
    #print("EC cleanup:", '%1.1f' % (cleanup_EC_end - cleanup_EC_start), "s")
    #print("-------------------------------------------------")
    print("Outputs:", '%1.1f' % (Cytoscape_end - Cytoscape_start - (cleanup_cytoscape_end - cleanup_cytoscape_start)), "s")
    print("Outputs cleanup:", '%1.1f' % (cleanup_cytoscape_end - cleanup_cytoscape_start), "s")
    

if __name__ == "__main__":
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline"
                                        "Version 1.0 © 2018")

    parser.add_argument("-c", "--config",   type=str,   help="Path to the configureation file")
    parser.add_argument("-1", "--pair1",    type=str,   help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2",    type=str,   help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single",   type=str,   help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, required=True, help="Path of the folder for the output of the pipeline")
    parser.add_argument("-t", "--num_threads", type=int, help="Maximum number of threads used by the pipeline")
    parser.add_argument("--nhost", action='store_true', help="Skip the host read removal step of the pipeline")
    parser.add_argument("--verbose_mode", type=str, help = "Decide how to handle the interim files, Compress them, or leave them alone.  Values are: keep, compress, quiet")
    parser.add_argument("-it", "--infernal_threads", type = int, help = "number of threads allowed for rRNA")
    parser.add_argument("-keep_second_split", action = 'store_true', help = "keep the interim second-split rRNA data")
    args = parser.parse_args()

    if (args.pair1 and not args.pair2) or (args.pair2 and not args.pair1):
        print("You must specify both forward and reverse reads for a paired-end run")
        sys.exit()
    elif args.single and (args.pair1 or args.pair2):
        print("You cannot specify both paired-end and single-end reads in a single run.")
        sys.exit()

    config_file     = args.config if args.config else ""
    pair_1          = args.pair1 if args.pair1 else ""
    pair_2          = args.pair2 if args.pair2 else ""
    single          = args.single if args.single else ""
    output_folder   = args.output_folder
    num_threads     = args.num_threads if args.num_threads else 0
    no_host         = args.nhost if args.nhost else False
    verbose_mode    = args.verbose_mode if args.verbose_mode else "quiet"
    keep_second_split = args.keep_second_split if args.keep_second_split else False
    infernal_limit  = args.infernal_threads if args.infernal_threads else 40
    if not (os.path.exists(output_folder)):
        print("output folder does not exist.  Now building directory.")
        os.makedirs(output_folder)
    os.chdir(output_folder)

    config = ConfigParser(interpolation = ExtendedInterpolation())
    if args.config:
        config.read(config_file)
        if not args.pair1 and not args.pair2 and not args.single:
            pair_1 = config["Sequences"]["pair1"] if config["Sequences"]["pair1"] else ""
            pair_2 = config["Sequences"]["pair2"] if config["Sequences"]["pair2"] else ""
            single = config["Sequences"]["single"] if config["Sequences"]["single"] else ""

    if pair_1 == "" and pair_2 == "" and single == "":
        print("You must specify paired-end or single-end reads as input for the pipeline.")
        sys.exit()

    args_pack = dict()
    args_pack["no_host"] = no_host
    args_pack["verbose_mode"] = verbose_mode
    args_pack["keep_second_split"] = keep_second_split
    args_pack["infernal_limit"] = infernal_limit
    
    main(config_file, pair_1, pair_2, single, output_folder, num_threads, args_pack)
