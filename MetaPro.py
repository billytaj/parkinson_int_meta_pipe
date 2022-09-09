#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.



#!/usr/bin/env python
from curses import meta
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import MetaPro_utilities as mpu
import MetaPro_stages as mps
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q


def main(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode):

    metapro_stage_obj = mps.mp_stage(config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, threads, args_pack, tutorial_mode)

    # This is the format we use to launch each stage of the pipeline.
    # We start a multiprocess that starts a subprocess.
    # The subprocess is created from the commands object

    # The quality filter stage
    #------------------------------------------------------------------------------
    metapro_stage_obj.mp_quality_filter()

    # The host read filter stage
    #-------------------------------------------------------------------------
    metapro_stage_obj.mp_host_filter()
        
    #-----------------------------------------------------------------
    # The vector contaminant filter stage
    metapro_stage_obj.mp_vector_filter()
    

    # ----------------------------------------------
    # rRNA removal stage
    metapro_stage_obj.mp_rRNA_filter()

    #---------------------------------------------------------------------------------------------------------------------
    # Duplicate repopulation
    metapro_stage_obj.mp_repop()
    # -------------------------------------------------------------
    # Assemble contigs
    metapro_stage_obj.mp_assemble()    

    # ----------------------------------------------
    # BWA gene annotation
    
    
    GA_BWA_start = time.time()
    GA_BWA_path = os.path.join(output_folder_path, GA_BWA_label)
    GA_BWA_jobs_folder = os.path.join(GA_BWA_path, "data", "jobs")
    #if not check_where_resume(GA_BWA_path, None, assemble_contigs_path):
    if mp_util.check_bypass_log(output_folder_path, GA_BWA_label):
        marker_path_list = []
        marker_file = "GA_split_fasta_contigs"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping", marker_file)
        else:
            job_name = "GA_prep_split_contigs"
            marker_path_list.append(marker_path)
            command_list = commands.create_split_ga_fasta_data_command(GA_BWA_label, assemble_contigs_label, "contigs", marker_file)
            mp_util.launch_and_create_with_mp_store(GA_BWA_label, job_name, commands, command_list)
        
        
        sections = ["singletons"]
        if(read_mode == "paired"):
            sections.extend(["pair_1", "pair_2"])
        for section in sections: 
            marker_file = "GA_split_fastq_" + section
            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping", marker_file)
            else:
                marker_path_list.append(marker_path)
                job_name = "GA_prep_split_" + section
                command_list = commands.create_split_ga_fastq_data_command(GA_BWA_label, assemble_contigs_label, section, marker_file)
                mp_util.launch_and_create_with_mp_store(GA_BWA_label, job_name, commands, command_list)
        mp_util.wait_for_mp_store()

        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_prep.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        
        #-------------------------------------------------------------------------
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                ref_path = paths.DNA_DB
                    
                command_list = ""
                if (ref_path.endswith(".fasta")):
                    ref_tag = os.path.basename(ref_path)
                    ref_tag = ref_tag.strip(".fasta")
                
                    file_tag = file_tag + "_" + ref_tag
                    job_name = "BWA" + "_" + file_tag
                    marker_file = file_tag + "_bwa"
                    marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                #this checker assumes that BWA only exports a file when it's finished running
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                    
                        #aug 10, 2021: new bigger chocophlan (from humann3) is in segments because we can't index it as a whole.  
                        #if the DB is still an old version, the tag should just say "chocophlan".  otherwise, it will say the chocophlan chunk name
                        
                        command_list = commands.create_BWA_annotate_command_v2(GA_BWA_label, ref_path, ref_tag, full_sample_path, marker_file)
                        mp_util.launch_and_create_with_hold(BWA_mem_threshold, BWA_job_limit, BWA_job_delay, GA_BWA_label, job_name, commands, command_list)
                else:
                    split_db = os.listdir(ref_path)
                    for db_segments in split_db:
                        if(db_segments.endswith(".fasta")):
                            segment_ref_path = os.path.join(ref_path, db_segments)
                            ref_tag = db_segments.strip(".fasta")
                            segment_file_tag = file_tag + "_" + ref_tag
                            job_name = "BWA" + "_" + segment_file_tag
                            marker_file = segment_file_tag + "_bwa"
                            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                            
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                                command_list = commands.create_BWA_annotate_command_v2(GA_BWA_label, segment_ref_path, ref_tag, full_sample_path, marker_file)
                                mp_util.launch_and_create_with_hold(BWA_mem_threshold, BWA_job_limit, BWA_job_delay, GA_BWA_label, job_name, commands, command_list)

        print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BWA_label)
            
    if mp_util.check_bypass_log(output_folder_path, GA_BWA_pp_label):
        marker_path_list = []
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                
                
                ref_path = paths.DNA_DB
                #chocophlan in many mutiple segments
                if (ref_path.endswith(".fasta")):
                    ref_tag = os.path.basename(ref_path)
                    ref_tag = ref_tag.strip(".fasta")
            
                
                    job_name = "BWA_pp" + "_" + file_tag + "_" + ref_tag
                    marker_file = file_tag + "_" + ref_tag +  "_bwa_pp"
                    marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                    
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = commands.create_BWA_pp_command_v2(GA_BWA_label, assemble_contigs_label, ref_tag, ref_path, full_sample_path, marker_file)
                        mp_util.launch_and_create_with_hold(BWA_pp_mem_threshold, BWA_pp_job_limit, BWA_pp_job_delay, GA_BWA_label, job_name, commands, command_list)
                        
                else:
                    #chocophlan in chunks
                    split_db = os.listdir(ref_path)
                    for db_segments in split_db:
                        if(db_segments.endswith(".fasta")):
                            segment_ref_path = os.path.join(ref_path, db_segments)
                            ref_tag = db_segments.strip(".fasta")
                            job_name = "BWA_pp" + "_" + file_tag + "_" + ref_tag
                            marker_file = file_tag + "_" + ref_tag + "_bwa_pp"
                            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                            
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                                command_list = commands.create_BWA_pp_command_v2(GA_BWA_label, assemble_contigs_label, ref_tag, segment_ref_path, full_sample_path, marker_file)
                                #print(dt.today(), "segmented BWA:", command_list)
                                #time.sleep(2)
                                mp_util.launch_and_create_with_hold(BWA_pp_mem_threshold, BWA_pp_job_limit, BWA_pp_job_delay, GA_BWA_label, job_name, commands, command_list)

                        
        print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
        mp_util.wait_for_mp_store()
        marker_file = "BWA_copy_contig_map"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:   
            marker_path_list.append(marker_path)
            command_list = commands.create_BWA_copy_contig_map_command(GA_BWA_label, assemble_contigs_label, marker_file)
            mp_util.launch_and_create_simple(GA_BWA_label, GA_BWA_label + "_copy_contig_map", commands, command_list)

        
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_pp.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BWA_pp_label)

    if mp_util.check_bypass_log(output_folder_path, GA_BWA_merge_label):
        #merge 
        marker_path_list = []
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
            
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                ref_path = paths.DNA_DB

                marker_file = file_tag + "_merge_fasta"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    job_name = "BWA_fasta_merge_" + file_tag
                    command_list = commands.create_merge_BWA_fasta_command(GA_BWA_label, full_sample_path, marker_file)
                    mp_util.launch_and_create_with_hold(BWA_pp_mem_threshold, BWA_pp_job_limit, BWA_pp_job_delay, GA_BWA_label, job_name, commands, command_list)

        print(dt.today(), "All BWA merge jobs have launched. waiting for sync")
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_merge.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BWA_merge_label)
        
 
    cleanup_GA_BWA_start = time.time()
    mp_util.delete_folder_simple(GA_BWA_jobs_folder)
    mp_util.clean_or_compress(GA_BWA_path, keep_all, keep_GA_BWA)
    
    cleanup_GA_BWA_end = time.time()
    GA_BWA_end = time.time()
    print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
    print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
    
    # ------------------------------------------------
    # BLAT gene annotation
    GA_BLAT_start = time.time()
    GA_BLAT_path = os.path.join(output_folder_path, GA_BLAT_label)
    GA_BLAT_data_folder = os.path.join(GA_BLAT_path, "data")
    GA_BLAT_jobs_folder = os.path.join(GA_BLAT_path, "data", "jobs")
    mp_util.make_folder(GA_BLAT_path)
    mp_util.make_folder(GA_BLAT_data_folder)
    mp_util.make_folder(GA_BLAT_jobs_folder)
    
    if mp_util.check_bypass_log(output_folder_path, GA_BLAT_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))

                delay_count = 0
                for fasta_db in os.listdir(paths.DNA_DB_Split):
                    if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                        job_name = "BLAT_" + file_tag + "_" + fasta_db
                        marker_file = file_tag + "_blat_" + fasta_db
                        marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                        blatout_path = os.path.join(GA_BLAT_path, "data", "0_blat", file_tag + "_"+fasta_db + ".blatout")
                        blat_queue_package = blatout_path+"|" + marker_file
                        
                        #This checker assume BLAT only exports a file when it's finished running
                        if(os.path.exists(marker_path)):
                            if(os.path.exists(blatout_path)):
                                #recover from a restart.  there will be files that have been missed.  thread would have deleted the file
                                print(dt.today(), "file still exists. adding to merge thread:", blatout_path)
                                #blat_file_queue.put(blat_queue_package)
                                
                            else:
                                print(dt.today(), "file doesn't exist anymore already merged", blatout_path)
                                
                            print(dt.today(), "BLAT job ran already, skipping:", marker_file)
                            #time.sleep(1)
                            continue
                            
                        else:
                            print(dt.today(), "RUNNING:", marker_file)
                            
                            marker_path_list.append(marker_path)
                            command_list = commands.create_BLAT_annotate_command_v2(GA_BLAT_label, full_sample_path, fasta_db, marker_file)
                            mp_util.launch_only_with_hold(BLAT_mem_threshold, BLAT_job_limit, BLAT_job_delay, job_name, commands, command_list)

        #---------------------------------------------------------------------------

        
        print(dt.today(), "final BLAT job removal. now waiting for mp-store flush")
        #note: this wait is disabled because we now have a separate thread.  it will hang if we enable it.
        print(dt.today(), "flushing mp_store")
        #mp_util.mp_store[:] = []        
        mp_util.wait_for_mp_store()
        print(dt.today(), "moving onto BLAT PP")
        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BLAT_label)
        
    #-------------------------------------------------
    #BLAT pp
    
    if mp_util.check_bypass_log(output_folder_path, GA_BLAT_pp_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                
                ref_path = paths.DNA_DB_Split  #the chocophlan chunks
                if (ref_path.endswith(".fasta")):
                    #single chocophlan mode
                    job_name = "BLAT_" + file_tag + "_pp"
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    marker_file = file_tag + "_blat_pp"
                    marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = commands.create_BLAT_pp_command_v2(GA_BLAT_label, full_sample_path, GA_BWA_label, ref_path, marker_file)
                        #mp_util.launch_and_create_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, GA_BLAT_label, job_name, commands, command_list)
                        mp_util.launch_only_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, job_name, commands, command_list)
                        
                else:
                    for fasta_db in os.listdir(ref_path):
                        if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                            #split chocophlan mode
                            #decode the chocophlan chunk, and supply the appropriate one.
                            #print("file tag:", file_tag.split("chocophlan"))
                            choco_chunk = fasta_db.split(".fasta")[0]
                            ref_file = os.path.join(ref_path, fasta_db)
                            #print("BLAT file tag:", file_tag, "|chunk: ", ref_file)
                            
                            job_name = "BLAT_" + file_tag + "_" + choco_chunk + "_pp"
                            full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                            #print("query file:", full_sample_path)
                            marker_file = file_tag + "_" + choco_chunk + "_blat_pp"
                            print("MARKER FILE:", marker_file)
                            marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                            
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                                command_list = commands.create_BLAT_pp_command_v3(GA_BLAT_label, full_sample_path, GA_BWA_label, ref_file, marker_file)
                                #print("Command list:", command_list)
                                #time.sleep(10)
                                #mp_util.launch_and_create_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, GA_BLAT_label, job_name, commands, command_list)
                                mp_util.launch_only_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, job_name, commands, command_list)
                            #time.sleep(10)

                
        print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
        mp_util.wait_for_mp_store()
        
        job_name = "GA_BLAT_copy_contigs"
        marker_file = "blat_copy_contig_map"
        marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_BLAT_copy_contig_map_command(GA_BLAT_label, GA_BWA_label, marker_file)
            mp_util.launch_and_create_simple(GA_BLAT_label, job_name, commands, command_list)
        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_pp.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BLAT_pp_label)
    #--------------------------------------------------------------
    # GA BLAT merge    
    if mp_util.check_bypass_log(output_folder_path, GA_BLAT_merge_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]

                marker_file = "BLAT_merge_" + file_tag
                marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping: ", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_BLAT_merge_fasta_command(GA_BLAT_label, file_tag, marker_file)
                    mp_util.launch_and_create_with_hold(BLAT_pp_mem_threshold, BLAT_pp_job_limit, BLAT_pp_job_delay, GA_BLAT_label, marker_file, commands, command_list)

        print(dt.today(), "submitted all BLAT merge jobs. waiting for sync")
        mp_util.wait_for_mp_store()

        final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_merge.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_BLAT_merge_label)

    #print(dt.today(), "stopping for a sanity check: BLAT merge")
    #sys.exit()

    cleanup_GA_BLAT_start = time.time()
    mp_util.delete_folder_simple(GA_BLAT_jobs_folder)
    mp_util.clean_or_compress(GA_BLAT_path, keep_all, keep_GA_BLAT)

    cleanup_GA_BLAT_end = time.time()
    GA_BLAT_end = time.time()
    print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
    print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
    
    # ------------------------------------------------------
    # Diamond gene annotation
    GA_DIAMOND_start = time.time()
    GA_DIAMOND_path = os.path.join(output_folder_path, GA_DIAMOND_label)
    GA_DIAMOND_tool_output_path = os.path.join(GA_DIAMOND_path, "data", "0_diamond")
    GA_DIAMOND_jobs_folder = os.path.join(GA_DIAMOND_path, "data", "jobs")
    #if not check_where_resume(None, GA_DIAMOND_tool_output_path, GA_BLAT_path, file_check_bypass = True):
    if mp_util.check_bypass_log(output_folder_path, GA_DIAMOND_label):
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_path)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_annotate_command_v2(GA_DIAMOND_label, full_sample_path, marker_file)
                    mp_util.launch_and_create_with_hold(DIAMOND_mem_threshold, DIAMOND_job_limit, DIAMOND_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                
        print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        mp_util.write_to_bypass_log(output_folder_path, GA_DIAMOND_label)
        
        
    #if not check_where_resume(GA_DIAMOND_path, None, GA_DIAMOND_tool_output_path, file_check_bypass = True):
    if mp_util.check_bypass_log(output_folder_path, GA_DIAMOND_pp_label):
        print(dt.today(), "DIAMOND PP threads used:", real_thread_count/2)
        marker_path_list = []
        for split_sample in os.listdir(os.path.join(GA_BLAT_path, "final_results")):
            if(split_sample.endswith(".fasta")):
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "DIAMOND_pp_" + file_tag
                full_sample_path = os.path.join(os.path.join(GA_BLAT_path, "final_results", split_sample))
                marker_file = file_tag + "_diamond_pp"
                marker_path = os.path.join(GA_DIAMOND_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = commands.create_DIAMOND_pp_command_v2(GA_DIAMOND_label, GA_BLAT_label, full_sample_path, marker_file)
                    mp_util.launch_and_create_with_hold(DIAMOND_pp_mem_threshold, DIAMOND_pp_job_limit, DIAMOND_pp_job_delay, GA_DIAMOND_label, job_name, commands, command_list)
                                    
        print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND_pp.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
        mp_util.write_to_bypass_log(output_folder_path, GA_DIAMOND_pp_label)
    
        
    
        cleanup_GA_DIAMOND_start = time.time()
        mp_util.delete_folder_simple(GA_DIAMOND_jobs_folder)
        mp_util.clean_or_compress(GA_DIAMOND_path, keep_all, keep_GA_DIAMOND)
        cleanup_GA_DIAMOND_end = time.time()
    GA_DIAMOND_end = time.time()
    print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
    print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
    
    
    GA_final_merge_start = time.time()
    GA_FINAL_MERGE_path = os.path.join(output_folder_path, GA_final_merge_label)
    if mp_util.check_bypass_log(output_folder_path, GA_final_merge_label):
        marker_file = "GA_final_merge"
        marker_path_p = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge_proteins")
        marker_path_m = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge_maps")
        marker_path_f = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge_fastq")
        if(os.path.exists(marker_path_p) and os.path.exists(marker_path_m) and os.path.exists(marker_path_f)):
            print(dt.today(), "skipping: GA final merge")
        else:
            command_list = commands.create_GA_final_merge_command(GA_final_merge_label, assemble_contigs_label, GA_BWA_label, GA_BLAT_label, GA_DIAMOND_label,  marker_file)
            job_name = "GA_final_merge"
            mp_util.subdivide_and_launch(GA_final_merge_job_delay, GA_final_merge_mem_threshold, GA_final_merge_job_limit, GA_final_merge_label, job_name, commands, command_list)
        
        #check if all_proteins.faa was generated
        all_proteins_path = os.path.join(output_folder_path, GA_final_merge_label, "final_results", "all_proteins.faa")
        if(os.path.exists(marker_path_p)):
            if(os.path.getsize(all_proteins_path) > 0):
                mp_util.write_to_bypass_log(output_folder_path, GA_final_merge_label)
                print(dt.today(), "All_proteins.faa is OK.  Continuing")
            else:
                sys.exit("GA final merge failed.  proteins weren't translated")
            
    GA_final_merge_end = time.time()
    print("GA final merge:", '%1.1f' % (GA_final_merge_end - GA_final_merge_start), "s")
    mp_util.clean_or_compress(GA_FINAL_MERGE_path, keep_all, keep_GA_final)
    
    # ------------------------------------------------------

    # Taxonomic annotation
    TA_start = time.time()
    TA_path = os.path.join(output_folder_path, taxon_annotation_label)
    TA_jobs_folder = os.path.join(TA_path, "data", "jobs")
    if mp_util.check_bypass_log(output_folder_path, taxon_annotation_label):
        #-----------------------------------------
        # stage 1
        marker_path_list = []
        #----------------------------------------------
        #centrifuge is too much of a RAM hog.  can't run more than 1 at a time
        sections = ["reads"]
        for section in sections:
            marker_file = "TA_centrifuge_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
                
        sections = ["contigs", "singletons"]
        if read_mode == "paired":
            sections.extend(["paired"])
            
        for section in sections:
            marker_file = "TA_kaiju_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_kaiju_command(taxon_annotation_label, assemble_contigs_label, section, marker_file)
                mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)        
        marker_file = "TA_taxon_pull"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_taxon_pull_command(taxon_annotation_label, GA_final_merge_label, marker_file)
            mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        print(dt.today(), "waiting for TA stage 1")
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(TA_path, "TA_stage_1.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        
        #--------------------------------------------------
        # stage 2
        marker_path_list = []
        sections = ["contigs"]
        for section in sections:
            marker_file = "TA_centrifuge_" + section
            marker_path = os.path.join(TA_jobs_folder, marker_file)
            
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:
                marker_path_list.append(marker_path)
                command_list = commands.create_TA_centrifuge_command(taxon_annotation_label, rRNA_filter_label, assemble_contigs_label, section, marker_file)
                mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        
        marker_file = "TA_kaiju_pp"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_file)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_kaiju_pp_command(taxon_annotation_label, marker_file)
            mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(TA_path, "TA_stage_2.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        #------------------------------------------------------------------

        #-----------------------------------------------------------------
        # stage 3
        marker_path_list = []
        marker_file = "TA_centrifuge_pp"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_centrifuge_pp_command(taxon_annotation_label, marker_file)
            mp_util.launch_and_create_with_hold(TA_mem_threshold, TA_job_limit, TA_job_delay, taxon_annotation_label, marker_file, commands, command_list)
        mp_util.wait_for_mp_store()
        final_checklist = os.path.join(TA_path, "TA_stage_3.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        #-----------------------------------------------
        # stage 4
        marker_path_list = []
        
        marker_file = "TA_final"
        marker_path = os.path.join(TA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            marker_path_list.append(marker_path)
            command_list = commands.create_TA_final_command(taxon_annotation_label, assemble_contigs_label, marker_file)
            mp_util.launch_and_create_simple(taxon_annotation_label, marker_file, commands, command_list)
        final_checklist = os.path.join(TA_path, "TA_final.txt")
        mp_util.check_all_job_markers(marker_path_list, final_checklist)
        
        if(os.path.exists(marker_path)):
            mp_util.write_to_bypass_log(output_folder_path, taxon_annotation_label)
            
    cleanup_TA_start = time.time()
    mp_util.clean_or_compress(TA_path, keep_all, keep_TA)
    cleanup_TA_end = time.time()
    TA_end = time.time()
    print("TA:", '%1.1f' % (TA_end - TA_start - (cleanup_TA_end - cleanup_TA_start)), "s")
    print("TA cleanup:", '%1.1f' % (cleanup_TA_end - cleanup_TA_start), "s")
    
    
    # ------------------------------------------------------
    # Detect EC annotation
    ec_annotation_path = os.path.join(output_folder_path, ec_annotation_label)
    EC_start = time.time()
    #There's a 2-step check.  We don't want it ti re-run either DETECT, or PRIAM+DIAMOND because they're too slow
    #if not check_where_resume(ec_annotation_path, None, GA_DIAMOND_path):
    #if check_bypass_log(output_folder_path, ec_annotation_label):
    EC_DETECT_start = time.time()
    ec_detect_path = os.path.join(ec_annotation_path, "data", "0_detect")
    #if not check_where_resume(job_label = None, full_path = ec_detect_path, dep_job_path = GA_DIAMOND_path):
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_detect_label):
        marker_file = "ec_detect"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            command_list = commands.create_EC_DETECT_command(ec_annotation_label, GA_final_merge_label, marker_file)
            mp_util.launch_and_create_with_mp_store(ec_annotation_label, marker_file, commands, command_list)
        
        
    EC_DETECT_end = time.time()
    print("EC DETECT:", '%1.1f' % (EC_DETECT_end - EC_DETECT_start), "s")
    
    # --------------------------------------------------------------
    # Priam EC annotation.  Why isn't it parallel? computing restraints.  Not enough mem
    EC_PRIAM_start = time.time()
    
    ec_priam_path = os.path.join(ec_annotation_path, "data", "1_priam")
    #if not check_where_resume(job_label = None, full_path = ec_priam_path, dep_job_path = GA_DIAMOND_path):
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_priam_label):
        marker_file = "ec_priam"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            command_list = commands.create_EC_PRIAM_command(ec_annotation_label, GA_final_merge_label, marker_file)
            if(os.path.exists(ec_priam_path)):
                print(dt.today(), "attempting PRIAM auto-resume")
                print("command:", command_list)
            mp_util.launch_and_create_with_mp_store(ec_annotation_label, marker_file, commands, command_list)
        
      
        #process.join()
    EC_PRIAM_end = time.time()
    print("EC PRIAM:", '%1.1f' % (EC_PRIAM_end - EC_PRIAM_start), "s")
    # --------------------------------------------------------------
    # DIAMOND EC annotation 
    EC_DIAMOND_start = time.time()
    ec_diamond_path = os.path.join(ec_annotation_path, "data", "2_diamond")
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        marker_file = "ec_diamond"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            job_name = "ec_diamond"
            command_list = commands.create_EC_DIAMOND_command(ec_annotation_label, GA_final_merge_label, marker_file)
            mp_util.launch_and_create_with_mp_store(ec_annotation_label, job_name, commands, command_list)
        
    EC_DIAMOND_end = time.time()
    mp_util.wait_for_mp_store()
    
    ec_detect_out   = os.path.join(ec_annotation_path, "data", "jobs", "ec_detect")
    ec_priam_out    = os.path.join(ec_annotation_path, "data", "jobs", "ec_priam")
    ec_diamond_out  = os.path.join(ec_annotation_path, "data", "jobs", "ec_diamond")
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_detect_label):
        if(os.path.exists(ec_detect_out)):
            mp_util.write_to_bypass_log(output_folder_path, ec_annotation_detect_label)
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_priam_label):
        if(os.path.exists(ec_priam_out)):
            mp_util.write_to_bypass_log(output_folder_path, ec_annotation_priam_label)
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_DIAMOND_label):
        if(os.path.exists(ec_diamond_out)):
            mp_util.write_to_bypass_log(output_folder_path, ec_annotation_DIAMOND_label)
    
    #----------------------------------------------------------------------
    # EC post process
    EC_post_start = time.time()
    #if not (check_where_resume(ec_annotation_path, None, GA_DIAMOND_path)):
    if mp_util.check_bypass_log(output_folder_path, ec_annotation_pp_label):
        
        marker_file = "ec_post"
        marker_path = os.path.join(ec_annotation_path, "data", "jobs", marker_file)
        
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping:", marker_file)
        else:
            command_list = commands.create_EC_postprocess_command(ec_annotation_label, GA_final_merge_label, marker_file)
            mp_util.launch_and_create_simple(ec_annotation_label, marker_file, commands, command_list)
        
        if(os.path.exists(marker_path)):
            mp_util.write_to_bypass_log(output_folder_path, ec_annotation_pp_label)
    
    cleanup_EC_start = time.time()
    mp_util.clean_or_compress(ec_annotation_path, keep_all, keep_EC)
    cleanup_EC_end = time.time()
    EC_post_end = time.time()
        
   
    EC_end = time.time()
    print("EC run:", '%1.1f' % (EC_end - EC_start), "s")
    print("EC cleanup:", '%1.1f' % (cleanup_EC_end - cleanup_EC_start), "s")
    
    # ------------------------------------------------------
    # RPKM Table and Cytoscape Network
    Cytoscape_start = time.time()
    network_path = os.path.join(output_folder_path, output_label)
    #if not check_where_resume(network_path, None, ec_annotation_path):
    
    if mp_util.check_bypass_log(output_folder, output_label):
        
        #phase 1
        if mp_util.check_bypass_log(output_folder, output_copy_gene_map_label):
            job_name = output_copy_gene_map_label
            command_list = commands.create_output_copy_gene_map_command(output_label, GA_final_merge_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
            
        if mp_util.check_bypass_log(output_folder, output_copy_taxa_label):
            job_name = output_copy_taxa_label
            command_list = commands.create_output_copy_taxa_command(output_label, taxon_annotation_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
        
        if mp_util.check_bypass_log(output_folder, output_contig_stats_label):
            job_name = output_contig_stats_label
            command_list = commands.create_output_contig_stats_command(output_label, assemble_contigs_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
            
        
            
        if not(no_host):
            print(dt.today(), "repopulating hosts for output")
            if mp_util.check_bypass_log(output_folder, output_unique_hosts_singletons_label):
                job_name = output_unique_hosts_singletons_label
                command_list = commands.create_output_unique_hosts_singletons_command(output_label, quality_filter_label, host_filter_label)
                mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
            
            if(read_mode == "paired"):
                if mp_util.check_bypass_log(output_folder, output_unique_hosts_pair_1_label):
                    job_name = output_unique_hosts_pair_1_label
                    command_list = commands.create_output_unique_hosts_pair_1_command(output_label, quality_filter_label, host_filter_label)
                    mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
                    
                if mp_util.check_bypass_log(output_folder, output_unique_hosts_pair_2_label):
                    job_name = output_unique_hosts_pair_2_label
                    command_list = commands.create_output_unique_hosts_pair_2_command(output_label, quality_filter_label, host_filter_label)
                    mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
                    
                    
        #repop vectors
        if mp_util.check_bypass_log(output_folder, output_unique_vectors_singletons_label):
            job_name = output_unique_vectors_singletons_label
            command_list = commands.create_output_unique_vectors_singletons_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
        
        if(read_mode == "paired"):
            if mp_util.check_bypass_log(output_folder, output_unique_vectors_pair_1_label):
                job_name = output_unique_vectors_pair_1_label
                command_list = commands.create_output_unique_vectors_pair_1_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
                mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
                
            if mp_util.check_bypass_log(output_folder, output_unique_vectors_pair_2_label):
                job_name = output_unique_vectors_pair_2_label
                command_list = commands.create_output_unique_vectors_pair_2_command(output_label, quality_filter_label, host_filter_label, vector_filter_label)
                mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
                
        print(dt.today(), "output report phase 1 launched.  waiting for sync")
        mp_util.wait_for_mp_store()
        
        mp_util.conditional_write_to_bypass_log(output_per_read_scores_label, "outputs/final_results", "input_per_seq_quality_report.csv")
        mp_util.conditional_write_to_bypass_log(output_copy_gene_map_label, "outputs/final_results", "final_gene_map.tsv")
        mp_util.conditional_write_to_bypass_log(output_copy_taxa_label, "outputs/final_results", "taxa_classifications.tsv")
        mp_util.conditional_write_to_bypass_log(output_contig_stats_label, "outputs/final_results", "contig_stats.txt")
        mp_util.conditional_write_to_bypass_log(output_unique_vectors_singletons_label, "outputs/data/4_full_vectors", "singletons_full_vectors.fastq")
        if(read_mode == "paired"):
            mp_util.conditional_write_to_bypass_log(output_unique_vectors_pair_1_label, "outputs/data/4_full_vectors", "pair_1_full_vectors.fastq")
            mp_util.conditional_write_to_bypass_log(output_unique_vectors_pair_2_label, "outputs/data/4_full_vectors", "pair_2_full_vectors.fastq")
            
        if not (no_host):
            mp_util.conditional_write_to_bypass_log(output_unique_hosts_singletons_label, "outputs/data/2_full_hosts", "singletons_full_hosts.fastq")
            if(read_mode == "paired"):
                mp_util.conditional_write_to_bypass_log(output_unique_hosts_pair_1_label, "outputs/data/2_full_hosts", "pair_1_full_hosts.fastq")
                mp_util.conditional_write_to_bypass_log(output_unique_hosts_pair_2_label, "outputs/data/2_full_hosts", "pair_2_full_hosts.fastq")
        #----------------------------------------------------------------------------
        #Phase 2
        if mp_util.check_bypass_log(output_folder, output_network_gen_label):
            command_list = commands.create_output_network_generation_command(output_label, GA_final_merge_label, taxon_annotation_label, ec_annotation_label)
            mp_util.launch_and_create_with_mp_store(output_label, output_network_gen_label, commands, command_list)
            
        if mp_util.check_bypass_log(output_folder, output_taxa_groupby_label):
            command_list = commands.create_output_taxa_groupby_command(output_label)
            mp_util.launch_and_create_with_mp_store(output_label, output_taxa_groupby_label, commands, command_list)
       
        print(dt.today(), "output report phase 2 launched.  waiting for sync")
        mp_util.wait_for_mp_store()
        mp_util.conditional_write_to_bypass_log(output_network_gen_label, "outputs/final_results", "RPKM_table.tsv")
        
        
        #-------------------------------------------------------------------
        #Phase 3
        if mp_util.check_bypass_log(output_folder, output_read_count_label):
            job_name = output_read_count_label
            command_list = commands.create_output_read_count_command(output_label, quality_filter_label, repop_job_label, GA_final_merge_label, ec_annotation_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
                               

        if mp_util.check_bypass_log(output_folder, output_per_read_scores_label):
            job_name = output_per_read_scores_label
            command_list = commands.create_output_per_read_scores_command(output_label, quality_filter_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)
            
        if mp_util.check_bypass_log(output_folder, output_ec_heatmap_label):
            job_name = output_ec_heatmap_label
            command_list = commands.create_output_EC_heatmap_command(output_label)
            mp_util.launch_and_create_with_mp_store(output_label, job_name, commands, command_list)    
        
        print(dt.today(), "output report phase 3 launched.  waiting for sync")
        mp_util.wait_for_mp_store()
        mp_util.conditional_write_to_bypass_log(output_read_count_label, "outputs/final_results", "read_count.tsv")
        mp_util.conditional_write_to_bypass_log(output_ec_heatmap_label, "outputs/final_results", "EC_coverage.csv")

        
    cleanup_cytoscape_start = time.time()
    mp_util.clean_or_compress(network_path, keep_all, keep_outputs)
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
    print("METAPRO metatranscriptomic analysis pipeline")
    # This is where the code starts
    # There's a few operating modes, mainly "docker", and "singularity".  These modes edit the pipeline filepaths

    parser = ArgumentParser(description="MetaPro - Meta-omic sequence processing and analysis pipeline"
                                        "Version 1.0.0 © 2022")

    parser.add_argument("-c", "--config",   type=str,   help="Path to the configureation file")
    parser.add_argument("-1", "--pair1",    type=str,   help="Path to the file containing the forward paired-end reads in fastq format")
    parser.add_argument("-2", "--pair2",    type=str,   help="Path to the file containing the reverse paired-end reads in fastq format")
    parser.add_argument("-s", "--single",   type=str,   help="Path to the file containing the single-end reads in fastq format")
    parser.add_argument("-con", "--contig",   type=str,   help="Tutorial use only: Path to the file containing the contig reads in fastq format")
    parser.add_argument("-o", "--output_folder", type=str, required=True, help="Path of the folder for the output of the pipeline")
    parser.add_argument("-t", "--num_threads", type=int, help="Maximum number of threads used by the pipeline")
    parser.add_argument("--nhost", "--no-host", action='store_true', help="Skip the host read removal step of the pipeline")
    parser.add_argument("--verbose_mode", type=str, help = "Decide how to handle the interim files, Compress them, or leave them alone.  Values are: keep, compress, quiet")
    parser.add_argument("--tutorial", type = str, help = "tutorial operating mode for MetaPro")
    
    args = parser.parse_args()
    
    config_file     = args.config if args.config else ""
    contig          = args.contig if args.contig else "None"
    pair_1          = args.pair1 if args.pair1 else ""
    pair_2          = args.pair2 if args.pair2 else ""
    single          = args.single if args.single else ""
    output_folder   = args.output_folder
    num_threads     = args.num_threads if args.num_threads else 0
    no_host         = args.nhost if args.nhost else False
    verbose_mode    = args.verbose_mode if args.verbose_mode else "quiet"
    tutorial_mode   = args.tutorial if args.tutorial else "none"

    if(tutorial_mode == "none"):
        if (args.pair1 and not args.pair2) or (args.pair2 and not args.pair1):
            print("You must specify both forward and reverse reads for a paired-end run")
            sys.exit()
        elif args.single and (args.pair1 or args.pair2):
            print("You cannot specify both paired-end and single-end reads in a single run.")
            sys.exit()

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
    
    print("=====================================")
    print("no-host:", no_host)
    print("verbose_mode:", verbose_mode)

    if (tutorial_mode != "none"):
        print("working in tutorial mode:", tutorial_mode)
        tutorial_main(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode)
    
    else:
        main(config_file, pair_1, pair_2, single, contig, output_folder, num_threads, args_pack, tutorial_mode)
