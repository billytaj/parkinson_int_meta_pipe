#!/usr/bin/env python
import sys
import os
import os.path
from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import MetaPro_commands as mpcom
import MetaPro_paths as mpp
import MetaPro_utilities as mpu
import MetaPro_tutorial as mpt
import time
import zipfile
import pandas as pd
import shutil
from datetime import datetime as dt
import psutil as psu
import threading as th
import queue as q

#stores code for stage-launch.
#makes for a neat package/capsule

class mp_stage:
    def __init__ (self, config_path, pair_1_path, pair_2_path, single_path, contig_path, output_folder_path, read_mode, threads, args_pack, tutorial_mode_string = None):
        #make our util obj
        #refresher: self -> instance var.  not self: class var (shared among class obj instances)
        tutorial_string = tutorial_mode_string
        self.mp_util = mpu.mp_util(self.output_folder_path)
        paths = mpp.tool_path_obj(config_path)
        self.segmented_chocophlan_flag = True
        if(paths.DNA_DB.endswith(".fasta")):
            self.segmented_chocophlan_flag = False
        self.no_host = args_pack["no_host"]
        self.verbose_mode = args_pack["verbose_mode"]
        self.rRNA_chunks = int(paths.rRNA_chunksize)
        self.GA_chunksize = int(paths.GA_chunksize)
        self.config_path = config_path
        self.pair_1_path = pair_1_path
        self.pair_2_path = pair_2_path
        self.single_path = single_path
        self.contig_path = contig_path
        self.read_mode = read_mode
        
        self.BWA_mem_threshold              = int(paths.BWA_mem_threshold)
        self.BLAT_mem_threshold             = int(paths.BLAT_mem_threshold)
        self.DIAMOND_mem_threshold          = int(paths.DIAMOND_mem_threshold)
        self.BWA_pp_mem_threshold           = int(paths.BWA_pp_mem_threshold)
        self.BLAT_pp_mem_threshold          = int(paths.BLAT_pp_mem_threshold)
        self.DIAMOND_pp_mem_threshold       = int(paths.DIAMOND_pp_mem_threshold)
        self.Infernal_mem_threshold         = int(paths.Infernal_mem_threshold)
        self.Barrnap_mem_threshold          = int(paths.Barrnap_mem_threshold)
        self.DETECT_mem_threshold           = int(paths.DETECT_mem_threshold)
        self.TA_mem_threshold               = int(paths.TA_mem_threshold)
        self.repop_mem_threshold            = int(paths.repop_mem_threshold)
        self.GA_final_merge_mem_threshold   = int(paths.GA_final_merge_mem_threshold)
        
        
        
        self.BWA_job_limit               = int(paths.BWA_job_limit)
        self.BLAT_job_limit              = int(paths.BLAT_job_limit)
        self.DIAMOND_job_limit           = int(paths.DIAMOND_job_limit)
        self.BWA_pp_job_limit            = int(paths.BWA_pp_job_limit)
        self.BLAT_pp_job_limit           = int(paths.BLAT_pp_job_limit)
        self.DIAMOND_pp_job_limit        = int(paths.DIAMOND_pp_job_limit)
        self.Infernal_job_limit          = int(paths.Infernal_job_limit)
        self.Barrnap_job_limit           = int(paths.Barrnap_job_limit)
        self.DETECT_job_limit            = int(paths.DETECT_job_limit)
        self.TA_job_limit                = int(paths.TA_job_limit)
        self.repop_job_limit             = int(paths.repop_job_limit)
        self.GA_final_merge_job_limit    = int(paths.GA_final_merge_job_limit)
        
        self.Infernal_job_delay          = float(paths.Infernal_job_delay)
        self.Barrnap_job_delay           = float(paths.Barrnap_job_delay)
        self.BWA_job_delay               = float(paths.BWA_job_delay)
        self.BLAT_job_delay              = float(paths.BLAT_job_delay)
        self.DIAMOND_job_delay           = float(paths.DIAMOND_job_delay)
        self.BWA_pp_job_delay            = float(paths.BWA_pp_job_delay)
        self.BLAT_pp_job_delay           = float(paths.BLAT_pp_job_delay)
        self.DIAMOND_pp_job_delay        = float(paths.DIAMOND_pp_job_delay)
        self.DETECT_job_delay            = float(paths.DETECT_job_delay)
        self.TA_job_delay                = float(paths.TA_job_delay)
        self.repop_job_delay             = float(paths.repop_job_delay)
        self.GA_final_merge_job_delay    = float(paths.GA_final_merge_job_delay)

        self.filter_stringency = paths.filter_stringency
        
        #-----------------------------------------------------
        self.keep_all                = paths.keep_all
        self.keep_quality            = paths.keep_quality
        self.keep_vector             = paths.keep_vector
        self.keep_host               = paths.keep_host
        self.keep_rRNA               = paths.keep_rRNA
        self.keep_repop              = paths.keep_repop
        self.keep_assemble_contigs   = paths.keep_assemble_contigs
        self.keep_GA_BWA             = paths.keep_GA_BWA
        self.keep_GA_BLAT            = paths.keep_GA_BLAT
        self.keep_GA_DIAMOND         = paths.keep_GA_DIAMOND
        self.keep_GA_final           = paths.keep_GA_final
        self.keep_TA                 = paths.keep_TA
        self.keep_EC                 = paths.keep_EC
        self.keep_outputs            = paths.keep_outputs
        
        #------------------------------------------------------------------------
        
        self.BWA_cigar_cutoff        = paths.BWA_cigar_cutoff
        self.BLAT_identity_cutoff    = paths.BLAT_identity_cutoff
        self.BLAT_length_cutoff      = paths.BLAT_length_cutoff
        self.BLAT_score_cutoff       = paths.BLAT_score_cutoff
        self.DIAMOND_identity_cutoff = paths.DIAMOND_identity_cutoff
        self.DIAMOND_length_cutoff   = paths.DIAMOND_length_cutoff
        self.DIAMOND_score_cutoff    = paths.DIAMOND_score_cutoff



        self.quality_filter_label                    = "quality_filter"
        self.host_filter_label                       = "host_read_filter"
        self.vector_filter_label                     = "vector_read_filter"
        self.rRNA_filter_label                       = "rRNA_filter"
        self.rRNA_filter_split_label                 = "rRNA_filter_split"
        self.rRNA_filter_convert_label               = "rRNA_filter_convert"
        self.rRNA_filter_barrnap_label               = "rRNA_filter_barrnap"
        self.rRNA_filter_barrnap_merge_label         = "rRNA_filter_barrnap_merge"
        self.rRNA_filter_barrnap_pp_label            = "rRNA_filter_barrnap_pp"
        self.rRNA_filter_infernal_label              = "rRNA_filter_infernal"
        self.rRNA_filter_infernal_prep_label         = "rRNA_filter_infernal_prep"
        self.rRNA_filter_splitter_label              = "rRNA_filter_splitter"
        self.rRNA_filter_post_label                  = "rRNA_filter_post"
        self.repop_job_label                         = "duplicate_repopulation"
        self.assemble_contigs_label                  = "assemble_contigs"
        self.destroy_contigs_label                   = "destroy_contigs"
        self.GA_BWA_label                            = "GA_BWA"
        self.GA_BWA_pp_label                         = "GA_BWA_pp"
        self.GA_BWA_merge_label                      = "GA_BWA_merge"
        self.GA_BLAT_label                           = "GA_BLAT"
        self.GA_BLAT_cleanup_label                   = "GA_BLAT_cleanup"
        self.GA_BLAT_cat_label                       = "GA_BLAT_cat"
        self.GA_BLAT_pp_label                        = "GA_BLAT_pp"
        self.GA_BLAT_merge_label                     = "GA_BLAT_merge"
        self.GA_DIAMOND_label                        = "GA_DIAMOND"
        self.GA_DIAMOND_pp_label                     = "GA_DIAMOND_pp"
        self.GA_final_merge_label                    = "GA_FINAL_MERGE"
        self.taxon_annotation_label                  = "taxonomic_annotation"
        self.ec_annotation_label                     = "enzyme_annotation"
        self.ec_annotation_detect_label              = "enzyme_annotation_detect"
        self.ec_annotation_priam_label               = "enzyme_annotation_priam"
        self.ec_annotation_DIAMOND_label             = "enzyme_annotation_DIAMOND"
        self.ec_annotation_pp_label                  = "enzyme_annotation_pp"
        self.output_label                            = "outputs"
        self.output_copy_gene_map_label              = "output_copy_gene_map"
        self.output_clean_EC_label                   = "output_clean_ec"
        self.output_copy_taxa_label                  = "output_copy_taxa"
        self.output_network_gen_label                = "output_network_generation"
        self.output_unique_hosts_singletons_label    = "output_unique_hosts_singletons"
        self.output_unique_hosts_pair_1_label        = "output_unique_hosts_pair_1"
        self.output_unique_hosts_pair_2_label        = "output_unique_hosts_pair_2"
        self.output_unique_vectors_singletons_label  = "output_unique_vectors_singletons"
        self.output_unique_vectors_pair_1_label      = "output_unique_vectors_pair_1"
        self.output_unique_vectors_pair_2_label      = "output_unique_vectors_pair_2"
        self.output_combine_hosts_label              = "output_combine_hosts"
        self.output_per_read_scores_label            = "output_per_read_scores"
        self.output_contig_stats_label               = "output_contig_stats"
        self.output_ec_heatmap_label                 = "output_ec_heatmap"
        self.output_taxa_groupby_label               = "output_taxa_groupby"
        self.output_read_count_label                 = "output_read_count"


        

        #timing vars
    
        self.quality_start                  = 0
        self.quality_end                    = 0
        self.cleanup_quality_start          = 0
        self.cleanup_quality_end            = 0
        
        self.host_start                     = 0
        self.host_end                       = 0
        self.cleanup_host_start             = 0
        self.cleanup_host_end               = 0
        
        self.vector_start                   = 0
        self.vector_end                     = 0
        self.cleanup_vector_start           = 0
        self.cleanup_vector_end             = 0
        
        self.rRNA_filter_start              = 0  
        self.rRNA_filter_end                = 0
        self.cleanup_rRNA_filter_start      = 0
        self.cleanup_rRNA_filter_end        = 0   
        
        self.repop_start                    = 0
        self.repop_end                      = 0
        self.cleanup_repop_start            = 0
        self.cleanup_repop_end              = 0
        
        self.assemble_contigs_start         = 0
        self.assemble_contigs_end           = 0
        self.cleanup_assemble_contigs_start = 0
        self.cleanup_assemble_contigs_end   = 0
        
        self.destroy_contigs_start          = 0
        self.destroy_contigs_end            = 0
        self.cleanup_destroy_contigs_start  = 0 
        self.cleanup_destroy_contigs_end    = 0
        
        self.GA_BWA_start                   = 0
        self.GA_BWA_end                     = 0
        self.cleanup_GA_BWA_start           = 0
        self.cleanup_GA_BWA_end             = 0
        
        self.GA_BLAT_start                  = 0
        self.GA_BLAT_end                    = 0
        self.cleanup_GA_BLAT_start          = 0
        self.cleanup_GA_BLAT_end            = 0
        
        self.GA_DIAMOND_start               = 0
        self.GA_DIAMOND_end                 = 0
        self.cleanup_GA_DIAMOND_start       = 0
        self.cleanup_GA_DIAMOND_end         = 0
        
        self.TA_start                       = 0
        self.TA_end                         = 0
        self.cleanup_TA_start               = 0
        self.cleanup_TA_end                 = 0
        
        self.EC_start                       = 0
        self.EC_end                         = 0

        self.EC_DETECT_start                = 0  
        self.EC_DETECT_end                  = 0
        
        self.EC_PRIAM_start                 = 0
        self.EC_PRIAM_end                   = 0
        
        self.EC_DIAMOND_start               = 0
        self.EC_DIAMOND_end                 = 0
        
        self.cleanup_EC_start               = 0
        self.cleanup_EC_end                 = 0
        
        self.Cytoscape_start                = 0
        self.Cytoscape_end                  = 0
        self.cleanup_cytoscape_start        = 0
        self.cleanup_cytoscape_end          = 0
        


        if not single_path == "":
            read_mode = "single"
            quality_encoding = mp_util.determine_encoding(single_path)
            print("ENCODING USED:", quality_encoding)
            print("OPERATING IN SINGLE-ENDED MODE")
        else:
            read_mode = "paired"
            quality_encoding = mp_util.determine_encoding(pair_1_path)
            print("ENCODING USED:", quality_encoding)
            print("OPERATING IN PAIRED-MODE")
            
        #number of threads to use/limit
        self.read_thread_count = threads
        if threads == 0:
            real_thread_count = mp.cpu_count()
        
            
        
        if(real_thread_count == 1):
            real_thread_count = 2
        print("number of threads used:", real_thread_count)         
                
        mp_store = []  # stores the multiprocessing processes

        # Creates our command object, for creating shellscripts.

        if self.read_mode == "single":
            self.commands = mpcom.mt_pipe_commands(self.no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = None, sequence_path_1=None, sequence_path_2=None, sequence_single=single_path, sequence_contigs = None)
        elif self.read_mode == "paired":
            self.commands = mpcom.mt_pipe_commands(self.no_host, Config_path=config_path, Quality_score=quality_encoding, Thread_count=real_thread_count, tutorial_keyword = None, sequence_path_1=pair_1_path, sequence_path_2=pair_2_path, sequence_single=None, sequence_contigs = None)
    

        #--------------------------------------------------------
        #working paths
        self.quality_path           = os.path.join(self.output_folder_path, self.quality_filter_label)
        self.host_path              = os.path.join(self.output_folder_path, self.host_filter_label)
        self.vector_path            = os.path.join(self.output_folder_path, self.vector_filter_label)
        self.rRNA_filter_path       = os.path.join(self.output_folder_path, self.rRNA_filter_label)
        self.repop_path             = os.path.join(self.output_folder_path, self.repop_job_label)
        self.assemble_contigs_path  = os.path.join(self.output_folder_path, self.assemble_contigs_label)
        
        

    #--------------------------------------------------------------------------------------------------------------
    # main calls
    def mp_quality_filter(self):
        self.quality_start = time.time()
        command_list = self.commands.create_quality_control_command(self.quality_filter_label)
        self.cleanup_quality_start, self.cleanup_quality_end = self.mp_util.launch_stage_simple(self.quality_filter_label, self.quality_path, self.commands, command_list, self.keep_all, self.keep_quality)
        self.quality_end = time.time()
        print("quality filter:", '%1.1f' % (self.quality_end - self.quality_start - (self.cleanup_quality_end - self.cleanup_quality_start)), "s")
        print("quality filter cleanup:", '%1.1f' %(self.cleanup_quality_end - self.cleanup_quality_start), "s")


    def mp_host_filter(self):
        if not self.no_host:
            self.host_start = time.time()
            #if not check_where_resume(host_path, None, quality_path):
            command_list = self.commands.create_host_filter_command(self.host_filter_label, self.quality_filter_label)
            self.cleanup_host_start, self.cleanup_host_end = self.mp_util.launch_stage_simple(self.host_filter_label, self.host_path, self.commands, command_list, self.keep_all, self.keep_host)
            self.host_end = time.time()
            print("host filter:", '%1.1f' % (self.host_end - self.host_start - (self.cleanup_host_end - self.cleanup_host_start)), "s")
            print("host filter cleanup:", '%1.1f' %(self.cleanup_host_end - self.cleanup_host_start),"s")

    def mp_vector_filter(self):
        self.vector_start = time.time()
        
        if self.no_host:
            #get dep args from quality filter
            #if not check_where_resume(vector_path, None, quality_path):
            command_list = self.commands.create_vector_filter_command(self.vector_filter_label, self.quality_filter_label)
            self.cleanup_vector_start, self.cleanup_vector_end = self.mp_util.launch_stage_simple(self.vector_filter_label, self.vector_path, self.commands, command_list, self.keep_all, self.keep_vector)

        else:
            #get the dep args from host filter
            #if not check_where_resume(vector_path, None, host_path):
            command_list = self.commands.create_vector_filter_command(self.vector_filter_label, self.host_filter_label)
            self.cleanup_vector_start, self.cleanup_vector_end = self.mp_util.launch_stage_simple(self.vector_filter_label, self.vector_path, self.commands, command_list, self.keep_all, self.keep_vector)
            
        self.vector_end = time.time()
        print("vector filter:", '%1.1f' % (self.vector_end - self.vector_start - (self.cleanup_vector_end - self.cleanup_vector_start)), "s")
        print("vector filter cleanup:", '%1.1f' % (self.cleanup_vector_end - self.cleanup_vector_start), "s")


    def mp_rRNA_filter(self):
        self.rRNA_filter_start = time.time()
        
        rRNA_filter_jobs_folder = os.path.join(self.rRNA_filter_path, "data", "jobs")
        #if not check_where_resume(self.rRNA_filter_path, None, vector_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_label): 
            marker_path_list = []
            sections = ["singletons"]
            if self.read_mode == "paired":
                sections.extend(["pair_1", "pair_2"])
            
            for section in reversed(sections):  #we go backwards due to a request by Ana.  pairs first, if applicable, then singletons
                #split the data, if necessary.
                #initial split -> by lines.  we can do both
                split_path = os.path.join(self.rRNA_filter_path, "data", section + "_fastq")
                barrnap_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section, section + "_barrnap")
                infernal_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section, section + "_infernal") 
                marker_file = "rRNA_filter_prep_" + section
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                #if not check_where_resume(job_label = None, full_path = second_split_path, dep_job_path = vector_path):
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_split_label + "_" + section):
                    print(dt.today(), "splitting:", section, " for rRNA filtration")
                    job_name = "rRNA_filter_prep_" + section
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_rRNA_filter_prep_command_v3(self.rRNA_filter_label, section, self.vector_filter_label, marker_file)
                    self.mp_util.launch_and_create_with_mp_store(self.rRNA_filter_label, job_name, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_prep.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            for element in sections:
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_split_label + "_" + element):
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_split_label + "_" + element)
                
            #-------------------------------------------------------------------------------------------------
            # Convert fastq segments to fasta
            
            for section in reversed(sections):
                split_path = os.path.join(self.rRNA_filter_path, "data", section + "_fastq")
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_convert_label + "_" + section):
                    marker_path_list = []
                    for item in os.listdir(split_path):
                        root_name = item.split(".")[0]
                        fasta_path = os.path.join(self.rRNA_filter_path, "data", section + "_fasta")
                        fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                        marker_file = root_name + "_convert_fasta"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        
                        fasta_out_size = os.stat(fasta_file).st_size if (os.path.exists(fasta_file)) else 0
                        if(fasta_out_size > 0) or (os.path.exists(marker_path)):
                            print(dt.today(), item, "already converted to fasta.  skipping")
                            continue
                        else:
                            job_name = root_name + "_convert_to_fasta"
                            marker_path_list.append(marker_path)
                            command_list = self.commands.create_rRNA_filter_convert_fastq_command("rRNA_filter", section, root_name+".fastq", marker_file)
                            self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                            
                    final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_convert_" + section + ".txt")
                    self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_convert_label + "_" + section)
                
                        
            #-------------------------------------------------------------------------------------------------
            # BARRNAP
            for section in reversed(sections):  
                #convert data to fasta, then run barrnap separately, then cat the barrnap, then run barrnap PP
                #split the data, if necessary.
                #initial split -> by lines.  we can do both
                split_path      = os.path.join(self.rRNA_filter_path, "data", section + "_fastq")
                fasta_path      = os.path.join(self.rRNA_filter_path, "data", section + "_fasta")
                barrnap_path    = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap")
                infernal_path   = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal") 
                
                mRNA_path       = os.path.join(self.rRNA_filter_path, "data", section + "_mRNA")
                
                #if not check_where_resume(job_label = None, full_path = barrnap_path, dep_job_path = vector_path):
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_barrnap_label + "_" + section):
                    concurrent_job_count = 0
                    batch_count = 0
                    marker_path_list = []
                    barrnap_org_list = ["arc", "bac", "euk", "mit"] #walk through all 4 types of organisms for barrnap
                    for item in os.listdir(fasta_path):
                        root_name = item.split(".")[0]
                        final_marker_file = root_name + "_barrnap_concat"
                        final_marker_path = os.path.join(rRNA_filter_jobs_folder, final_marker_file)
                        final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                        for barrnap_org in barrnap_org_list:
                            barrnap_out_file = os.path.join(barrnap_path, root_name + "_" + barrnap_org + ".barrnap_out")
                            fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                            fastq_file = os.path.join(split_path, root_name + ".fastq")
                            barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                            marker_file = root_name + "_barrnap_" + barrnap_org
                            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                            barrnap_out_size = os.stat(barrnap_out_file).st_size if (os.path.exists(barrnap_out_file)) else 0
                            
                            if(os.path.exists(final_marker_path)):
                                print(dt.today(), "skipping barrnap.  data already merged", final_marker_path)
                                continue
                            else:
                                if((barrnap_out_size > 0) and (os.path.exists(marker_path))):
                                    continue
                                else:
                                    marker_path_list.append(marker_path)
                                    command_list = ""
                                    if(barrnap_org == "arc"):
                                        command_list = self.commands.create_rRNA_filter_barrnap_arc_command("rRNA_filter", section, root_name, marker_file)
                                    elif(barrnap_org == "bac"):
                                        command_list = self.commands.create_rRNA_filter_barrnap_bac_command("rRNA_filter", section, root_name, marker_file)
                                    elif(barrnap_org == "euk"):
                                        command_list = self.commands.create_rRNA_filter_barrnap_euk_command("rRNA_filter", section, root_name, marker_file)
                                    elif(barrnap_org == "mit"):
                                        command_list = self.commands.create_rRNA_filter_barrnap_mit_command("rRNA_filter", section, root_name, marker_file)
                                    job_name = marker_file
                                    
                                    self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                                    
                    print(dt.today(), "waiting for Barrnap jobs to finish")
                    self.mp_util.wait_for_mp_store()

                    final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_" + section + ".txt")
                    self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                    
                    #------------------------------------------------------
                    #merge the barrnap data
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_barrnap_merge_label + "_" + section):
                        marker_path_list = []
                        for item in os.listdir(fasta_path):
                            root_name = item.split(".")[0]
                            marker_file = root_name + "_barrnap_cat"
                            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                            final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                            #final_barrnap_out_size  = os.stat(final_barrnap_out).st_size if (os.path.exists(final_barrnap_out)) else 0
                            
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "barrnap already merged. skipping:", item)
                                continue
                            else:
                                job_name = root_name + "_barrnap_cat"
                                marker_path_list.append(marker_path)
                                command_list = self.commands.create_rRNA_filter_barrnap_cat_command("rRNA_filter", section, root_name, marker_file)
                                self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                    print(dt.today(), "waiting for Barrnap pp to finish")
                    self.mp_util.wait_for_mp_store()
                    final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_cat_" + section  + ".txt")
                    self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_barrnap_merge_label + "_" + section)

                    #-----------------------------------------------------
                    #run the barrnap PP
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_barrnap_pp_label + "_" + section):
                        marker_path_list = []
                        for item in os.listdir(fasta_path):
                            root_name = item.split(".")[0]
                            barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                            barrnap_mRNA_out_size   = os.stat(barrnap_mrna_file).st_size if (os.path.exists(barrnap_mrna_file)) else 0
                            marker_file = root_name + "_barrnap_pp"
                            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "barrnap pp already run.  skipping:", item)
                                continue
                            else:
                                job_name = root_name + "_barrnap_pp"
                                marker_path_list.append(marker_path)
                                command_list = self.commands.create_rRNA_filter_barrnap_pp_command("rRNA_filter", section, root_name + ".fastq", marker_file)
                                self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                        
                        print(dt.today(), "waiting for Barrnap pp to finish")
                        self.mp_util.wait_for_mp_store()
                        final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_" + section +  ".txt")
                        self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                        self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_barrnap_label + "_" + section)
                
            #----------------------------------------------------------------------------
            # INFERNAL
            for section in reversed(sections):  
                #split the data, if necessary.
                #initial split -> by lines.  we can do both
                barrnap_mRNA_fastq_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap_mRNA")
                infernal_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal") 
                barrnap_mRNA_fasta_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap_mRNA_fasta")
                splitter_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal_mRNA")
            
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_infernal_prep_label + "_" + section):
                    concurrent_job_count = 0
                    batch_count = 0
                    #these jobs now have to be launched in segments
                    for item in os.listdir(barrnap_mRNA_fastq_path):
                    
                        if(item.endswith("_barrnap_mRNA.fastq")):
                            root_name = item.split(".")[0]
                            marker_file = root_name + "_infernal_prep"
                            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                            infernal_prep_out_file = os.path.join(barrnap_mRNA_fasta_path, root_name + ".fasta")
                            infernal_prep_file_size = os.stat(infernal_prep_out_file).st_size if (os.path.exists(infernal_prep_out_file)) else 0
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "Infernal prep already ran on this sample.  skipping", item)
                                continue
                            
                            else:
                                marker_path_list.append(marker_path)
                                job_name = "rRNA_filter_infernal_prep_" + root_name
                                command_list = self.commands.create_rRNA_filter_infernal_prep_command("rRNA_filter", section, item, root_name, marker_file)
                                self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
                                
                    print(dt.today(), "final batch: infernal prep")
                    self.mp_util.wait_for_mp_store()
                    final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_prep_" + section + ".txt")
                    self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_infernal_prep_label + "_" + section)
                

                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_infernal_label + "_" + section):
                    marker_path_list = []
                    for item in os.listdir(barrnap_mRNA_fasta_path):
                        #using a job marker is ineffective.  The marker will still write 
                        root_name = item.split("_barrnap_mRNA")[0]
                        marker_file = root_name + "_infernal"
                        marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                        
                        if(os.path.exists(marker_path)):
                            print(dt.today(), "infernal already run. skipping:", root_name + "_infernal")
                            continue
                        else:
                            marker_path_list.append(marker_path)
                            inf_command = self.commands.create_rRNA_filter_infernal_command("rRNA_filter", section, root_name, marker_file)
                            job_name = "rRNA_filter_infernal_" + root_name
                            #launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, self.commands, inf_command)
                            self.mp_util.launch_and_create_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, self.rRNA_filter_label, job_name, self.commands, inf_command)
                            
                            
                    print(dt.today(), "final batch: infernal")
                    self.mp_util.wait_for_mp_store()
                    final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_" + section + ".txt")
                    self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                    self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_infernal_label + "_" + section)
                
                if (section != "pair_2"):
                    if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_splitter_label + "_" + section):
                        marker_path_list = []
                        for item in os.listdir(barrnap_mRNA_fasta_path):
                            root_name = item.split("_barrnap_mRNA")[0]
                            splitter_out_file = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal_mRNA", root_name + "_mRNA.fastq")
                            splitter_out_file_size = os.stat(splitter_out_file).st_size if os.path.exists(splitter_out_file) else 0
                            marker_file = root_name + "_infernal_pp"
                            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "infernal mRNA splitter already run. skipping:", marker_file)
                                print("file size:", splitter_out_file_size, "file:", splitter_out_file)
                                continue
                            else:
                                job_name = "rRNA_filter_infernal_splitter_" + root_name
                                marker_path_list.append(marker_path)
                                command_list = self.commands.create_rRNA_filter_splitter_command("rRNA_filter", section, root_name, marker_file)
                                print(command_list)
                                self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
                                
                        print(dt.today(), "final batch: infernal splitter")
                        self.mp_util.wait_for_mp_store()
                        final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_splitter_" + section + ".txt")
                        self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                        self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_splitter_label + "_" + section)
                else:
                    print(dt.today(), "not calling Infernal rRNA splitter on pair 2.  data handled by pair 1 as a combination")
            
                        
                        
            marker_path_list = []
            for section in reversed(sections):
                if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_post_label + "_" + section):
                    print(dt.today(), "now running rRNA filter post:", section)
                    marker_file = section + "_rRNA_packup"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    job_name = "rRNA_post_cat"
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_rRNA_filter_final_cat_command("rRNA_filter", section, marker_file), self.assemble_contigs_labe
                    print("command list:", command_list)
                    self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
                    
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_final_cat.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            for section in reversed(sections):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_splitter_label + "_" + section)
            
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.rRNA_filter_label)
            self.cleanup_rRNA_filter_start = time.time()
            self.mp_util.delete_folder_simple(rRNA_filter_jobs_folder)
            self.mp_util.clean_or_compress(self.rRNA_filter_path, self.keep_all, self.keep_rRNA)
            
            self.cleanup_rRNA_filter_end = time.time()
        self.rRNA_filter_end = time.time()
        
        print("rRNA filter:", '%1.1f' % (self.rRNA_filter_end - self.rRNA_filter_start - (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start)), "s")
        print("rRNA filter cleanup:", '%1.1f' % (self.cleanup_rRNA_filter_end - self.cleanup_rRNA_filter_start), "s")

    def mp_repop(self):
        
        self.repop_start = time.time()
        #if not check_where_resume(repop_job_path, None, rRNA_filter_path):
        if self.mp_util.check_bypass_log(self.output_folder_path, self.repop_job_label):
            job_name = self.repop_job_label
            command_list = self.commands.create_repop_command_v2_step_1(self.repop_job_label, self.quality_filter_label, self.rRNA_filter_label)
            self.mp_util.subdivide_and_launch(self.repop_job_delay, self.repop_mem_threshold, self.repop_job_limit, self.repop_job_label, job_name, self.commands, command_list)
            self.mp_util.wait_for_mp_store()
            
            if(self.read_mode == "paired"):
                job_name = self.repop_job_label
                command_list = self.commands.create_repop_command_v2_step_2(self.repop_job_label, self.quality_filter_label, self.rRNA_filter_label)
                self.mp_util.subdivide_and_launch(self.repop_job_delay, self.repop_mem_threshold, self.repop_job_limit, self.repop_job_label, job_name, self.commands, command_list)
                self.mp_util.wait_for_mp_store()
            
            self.mp_util.write_to_bypass_log(self.output_folder_path, self.repop_job_label)
            
            self.cleanup_repop_start = time.time()
            self.mp_util.clean_or_compress(self.repop_path, self.keep_all, self.keep_repop)
            self.cleanup_repop_end = time.time()
            
        self.repop_end = time.time()
        print("repop:", '%1.1f' % (self.repop_end - self.repop_start - (self.cleanup_repop_end - self.cleanup_repop_start)), "s")
        print("repop cleanup:", '%1.1f' % (self.cleanup_repop_end - self.cleanup_repop_start), "s")

    def mp_assemble():
        self.assemble_contigs_start = time.time()
        
        
        #if not check_where_resume(assemble_contigs_path, None, repop_job_path):
        
        if self.mp_util.check_bypass_log(self.output_folder_path, self.assemble_contigs_label):
            job_name = self.assemble_contigs_label
            command_list = self.commands.create_assemble_contigs_command(self.assemble_contigs_label, self.repop_job_label)
            self.mp_util.launch_and_create_simple(self.assemble_contigs_label, job_name, self.commands, command_list)
            mgm_file = os.path.join(self.assemble_contigs_path, "data", "1_mgm", "gene_report.txt")
            if(os.path.exists(mgm_file)):
                self.mp_util.write_to_bypass_log(self.output_folder_path, self.assemble_contigs_label)
            else:
                sys.exit("mgm did not run.  look into it.  pipeline stopping here")
            
            self.cleanup_assemble_contigs_start = time.time()
            self.mp_util.clean_or_compress(self.assemble_contigs_path, self.keep_all, self.keep_assemble_contigs)

            self.cleanup_assemble_contigs_end = time.time()
        
        
        self.assemble_contigs_end = time.time()
        print("assemble contigs:", '%1.1f' % (self.assemble_contigs_end - self.assemble_contigs_start - (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start)), "s")    
        print("assemble contigs cleanup:", '%1.1f' % (self.cleanup_assemble_contigs_end - self.cleanup_assemble_contigs_start), "s")
    
    
    def mp_GA_BWA(self):

    #-------------------------------------------------------------------------------------------------------------
    #tutorial calls
    def tutorial_quality(self):
        print(dt.today(), "working on:", self.tutorial_string)
        command_list = self.commands.create_quality_control_command(self.quality_filter_label)
        job_name = self.quality_filter_label
        self.mp_util.launch_and_create_simple(self.quality_filter_label, job_name, self.commands, command_list)

    def tutorial_host(self):
        print(dt.today(), "working on:", self.tutorial_string)
        job_name = self.host_filter_label
        command_list = self.commands.create_host_filter_command(self.host_filter_label, self.quality_filter_label)
        self.mp_util.launch_and_create_simple(self.host_filter_label, job_name, self.commands, command_list)

    def tutorial_vectors(self):
        print(dt.today(), "working on:", self.tutorial_string)
        job_name = self.vector_filter_label
        command_list = self.commands.create_vector_filter_command(self.vector_filter_label, self.quality_filter_label)
        self.mp_util.launch_and_create_simple(self.vector_filter_label, job_name, self.commands, command_list)


    def tutorial_rRNA(self):
        print(dt.today(), "working on:", self.tutorial_string)
        
        rRNA_filter_jobs_folder = os.path.join(self.rRNA_filter_path, "data", "jobs")
        marker_path_list = []
        sections = ["singletons"]
        if self.read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        job_name = "rRNA_filter_prep_tutorial"
        command_list = self.commands.create_rRNA_filter_prep_command_v3(self.rRNA_filter_label, "tutorial", self.vector_filter_label, "tutorial")
        self.mp_util.launch_and_create_simple(self.rRNA_filter_label, job_name, self.commands, command_list)

        for section in reversed(sections):
            split_path = os.path.join(self.rRNA_filter_path, "data", section + "_fastq")
            if self.mp_util.check_bypass_log(self.output_folder_path, self.rRNA_filter_convert_label + "_" + section):
                marker_path_list = []
                for item in os.listdir(split_path):
                    root_name = item.split(".")[0]
                    fasta_path = os.path.join(self.rRNA_filter_path, "data", section + "_fasta")
                    fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                    marker_file = root_name + "_convert_fasta"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    
                    fasta_out_size = os.stat(fasta_file).st_size if (os.path.exists(fasta_file)) else 0
                    if(os.path.exists(marker_path)):
                        print(dt.today(), item, "already converted to fasta.  skipping")
                        continue
                    else:
                        job_name = root_name + "_convert_to_fasta"
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_rRNA_filter_convert_fastq_command("rRNA_filter", section, root_name+".fastq", marker_file)
                        self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                        
                        
        # BARRNAP
        for section in reversed(sections):  
            #convert data to fasta, then run barrnap separately, then cat the barrnap, then run barrnap PP
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            split_path      = os.path.join(self.rRNA_filter_path, "data", section + "_fastq")
            fasta_path      = os.path.join(self.rRNA_filter_path, "data", section + "_fasta")
            barrnap_path    = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap")
            infernal_path   = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal") 
            
            mRNA_path       = os.path.join(self.rRNA_filter_path, "data", section + "_mRNA")
            
            concurrent_job_count = 0
            batch_count = 0
            marker_path_list = []
            
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                final_marker_file = root_name + "_barrnap_cat"
                final_marker_path = os.path.join(rRNA_filter_jobs_folder, final_marker_file)
                barrnap_arc_out_file = os.path.join(barrnap_path, root_name + "_arc.barrnap_out")
                barrnap_bac_out_file = os.path.join(barrnap_path, root_name + "_bac.barrnap_out")
                barrnap_euk_out_file = os.path.join(barrnap_path, root_name + "_euk.barrnap_out")
                barrnap_mit_out_file = os.path.join(barrnap_path, root_name + "_mit.barrnap_out")
                final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                fasta_file = os.path.join(fasta_path, root_name + ".fasta")
                fastq_file = os.path.join(split_path, root_name + ".fastq")
                barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                marker_file_arc = root_name + "_barrnap_arc"
                marker_file_bac = root_name + "_barrnap_bac"
                marker_file_euk = root_name + "_barrnap_euk"
                marker_file_mit = root_name + "_barrnap_mit"
                
                marker_path_arc = os.path.join(rRNA_filter_jobs_folder, marker_file_arc)
                marker_path_bac = os.path.join(rRNA_filter_jobs_folder, marker_file_bac)
                marker_path_euk = os.path.join(rRNA_filter_jobs_folder, marker_file_euk)
                marker_path_mit = os.path.join(rRNA_filter_jobs_folder, marker_file_mit)
                
                barrnap_arc_out_size    = os.stat(barrnap_arc_out_file).st_size if (os.path.exists(barrnap_arc_out_file)) else 0
                barrnap_bac_out_size    = os.stat(barrnap_bac_out_file).st_size if (os.path.exists(barrnap_bac_out_file)) else 0
                barrnap_euk_out_size    = os.stat(barrnap_euk_out_file).st_size if (os.path.exists(barrnap_euk_out_file)) else 0
                barrnap_mit_out_size    = os.stat(barrnap_mit_out_file).st_size if (os.path.exists(barrnap_mit_out_file)) else 0
                
                if(os.path.exists(final_marker_path)):
                    print(dt.today(), "Job already run. skipping:", final_marker_file, final_marker_path)
                    time.sleep(5)
                    continue
                else:
                    if((barrnap_arc_out_size > 0) and (os.path.exists(marker_path_arc))):
                        print(dt.today(), "barrnap arc already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_arc"
                        marker_path_list.append(marker_path_arc)
                        command_list = self.commands.create_rRNA_filter_barrnap_arc_command("rRNA_filter", section, root_name, marker_file_arc)
                        self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                        
                        
                    if((barrnap_bac_out_size > 0) and (os.path.exists(marker_path_bac))):
                        print(dt.today(), "barrnap bac already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_bac"
                        marker_path_list.append(marker_path_bac)
                        command_list = self.commands.create_rRNA_filter_barrnap_bac_command("rRNA_filter", section, root_name, marker_file_bac)
                        self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                        
                    if((barrnap_euk_out_size > 0) and (os.path.join(marker_path_euk))):
                        print(dt.today(), "barrnap euk already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_euk"
                        marker_path_list.append(marker_path_euk)
                        command_list = self.commands.create_rRNA_filter_barrnap_euk_command("rRNA_filter", section, root_name, marker_file_euk)
                        self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
                        
                    if((barrnap_mit_out_size > 0) and (os.path.join(marker_path_mit))):
                        print(dt.today(), "barrnap mit already run.  skipping:", item) 
                        continue
                    else:
                        job_name = root_name + "_barrnap_mit"
                        marker_path_list.append(marker_path_mit)
                        command_list = self.commands.create_rRNA_filter_barrnap_mit_command("rRNA_filter", section, root_name, marker_file_mit)
                        self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
            print(dt.today(), "waiting for Barrnap jobs to finish")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_" + section + ".txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
        
            #------------------------------------------------------
            #merge the barrnap data
            marker_path_list = []
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                marker_file = root_name + "_barrnap_cat"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                final_barrnap_out    = os.path.join(barrnap_path, root_name + ".barrnap_out")
                final_barrnap_out_size  = os.stat(final_barrnap_out).st_size if (os.path.exists(final_barrnap_out)) else 0
                
                if((final_barrnap_out_size > 0) and (os.path.exists(marker_path))):
                    print(dt.today(), "barrnap already merged. skipping:", item)
                    continue
                else:
                    job_name = root_name + "_barrnap_cat"
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_rRNA_filter_barrnap_cat_command("rRNA_filter", section, root_name, marker_file)
                    self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
            print(dt.today(), "waiting for Barrnap pp to finish")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_cat_" + section  + ".txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
            #-----------------------------------------------------
            #run the barrnap PP
            marker_path_list = []
            for item in os.listdir(fasta_path):
                root_name = item.split(".")[0]
                barrnap_mrna_file   = os.path.join(mRNA_path, root_name + "_barrnap_mRNA.fastq")
                barrnap_mRNA_out_size   = os.stat(barrnap_mrna_file).st_size if (os.path.exists(barrnap_mrna_file)) else 0
                marker_file = root_name + "_barrnap_pp"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                if(barrnap_mRNA_out_size > 0):
                    print(dt.today(), "barrnap pp already run.  skipping:", item)
                    continue
                else:
                    job_name = root_name + "_barrnap_pp"
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_rRNA_filter_barrnap_pp_command("rRNA_filter", section, root_name + ".fastq", marker_file)
                    self.mp_util.launch_only_with_hold(self.Barrnap_mem_threshold, self.Barrnap_job_limit, self.Barrnap_job_delay, job_name, self.commands, command_list)
            
            print(dt.today(), "waiting for Barrnap pp to finish")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_barrnap_" + section +  ".txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                
        #----------------------------------------------------------------------------
        # INFERNAL
        for section in reversed(sections):  
            #split the data, if necessary.
            #initial split -> by lines.  we can do both
            barrnap_mRNA_fastq_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap_mRNA")
            infernal_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal") 
            barrnap_mRNA_fasta_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_barrnap_mRNA_fasta")
            splitter_path = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal_mRNA")
        
            concurrent_job_count = 0
            batch_count = 0
            #these jobs now have to be launched in segments
            for item in os.listdir(barrnap_mRNA_fastq_path):
               
                if(item.endswith("_barrnap_mRNA.fastq")):
                    root_name = item.split(".")[0]
                    marker_file = root_name + "_infernal_prep"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    infernal_prep_out_file = os.path.join(barrnap_mRNA_fasta_path, root_name + ".fasta")
                    infernal_prep_file_size = os.stat(infernal_prep_out_file).st_size if (os.path.exists(infernal_prep_out_file)) else 0
                    if((infernal_prep_file_size > 0) and (os.path.exists(marker_path))):
                        print(dt.today(), "Infernal prep already ran on this sample.  skipping", item)
                        continue
                    
                    else:
                        marker_path_list.append(marker_path)
                        job_name = "rRNA_filter_infernal_prep_" + root_name
                        command_list = self.commands.create_rRNA_filter_infernal_prep_command("rRNA_filter", section, item, root_name, marker_file)
                        self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
                        
            print(dt.today(), "final batch: infernal prep")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_prep_" + section + ".txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            

            marker_path_list = []
            for item in os.listdir(barrnap_mRNA_fasta_path):
                #using a job marker is ineffective.  The marker will still write 
                root_name = item.split("_barrnap_mRNA")[0]
                marker_file = root_name + "_infernal"
                marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                
                if(os.path.exists(marker_path)):
                    print(dt.today(), "infernal already run. skipping:", root_name + "_infernal")
                    continue
                else:
                    marker_path_list.append(marker_path)
                    inf_command = self.commands.create_rRNA_filter_infernal_command("rRNA_filter", section, root_name, marker_file)
                    job_name = "rRNA_filter_infernal_" + root_name
                    #launch_only_with_hold(mp_store, Infernal_mem_threshold, Infernal_job_limit, Infernal_job_delay, job_name, self.commands, inf_command)
                    self.mp_util.launch_and_create_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, self.rRNA_filter_label, job_name, self.commands, inf_command)
                    
                    
            print(dt.today(), "final batch: infernal")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_" + section + ".txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
            if (section != "pair_2"):
                marker_path_list = []
                for item in os.listdir(barrnap_mRNA_fasta_path):
                    root_name = item.split("_barrnap_mRNA")[0]
                    splitter_out_file = os.path.join(self.output_folder_path, self.rRNA_filter_label, "data", section + "_infernal_mRNA", root_name + "_mRNA.fastq")
                    splitter_out_file_size = os.stat(splitter_out_file).st_size if os.path.exists(splitter_out_file) else 0
                    marker_file = root_name + "_infernal_pp"
                    marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
                    if((splitter_out_file_size > 0) and (os.path.exists(marker_path))):
                        print(dt.today(), "infernal mRNA splitter already run. skipping:", marker_file)
                        print("file size:", splitter_out_file_size, "file:", splitter_out_file)
                        continue
                    else:
                        job_name = "rRNA_filter_infernal_splitter_" + root_name
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_rRNA_filter_splitter_command("rRNA_filter", section, root_name, marker_file)
                        print(command_list)
                        self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
                        
                print(dt.today(), "final batch: infernal splitter")
                self.mp_util.wait_for_mp_store()
                final_checklist = os.path.join(self.rRNA_filter_path, "rRNA_filter_infernal_splitter_" + section + ".txt")
                self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                
            else:
                print(dt.today(), "not calling Infernal rRNA splitter on pair 2.  data handled by pair 1 as a combination")
        
                    
                    
        marker_path_list = []
        for section in reversed(sections):
            print(dt.today(), "now running rRNA filter post:", section)
            marker_file = section + "_rRNA_packup"
            marker_path = os.path.join(rRNA_filter_jobs_folder, marker_file)
            job_name = "rRNA_post_cat"
            marker_path_list.append(marker_path)
            command_list = self.commands.create_rRNA_filter_final_cat_command("rRNA_filter", section, marker_file)
            print("command list:", command_list)
            self.mp_util.launch_only_with_hold(self.Infernal_mem_threshold, self.Infernal_job_limit, self.Infernal_job_delay, job_name, self.commands, command_list)
            
        self.mp_util.wait_for_mp_store

    def tutorial_repop(self):
        print(dt.today(), "working on:", self.tutorial_string)
        quality_path = os.path.join(self.output_folder_path, self.quality_filter_label)
        if(os.path.exists(quality_path)):
            job_name = self.repop_job_label
            command_list = self.commands.create_repop_command(self.repop_job_label, self.quality_filter_label, self.rRNA_filter_label)
            self.mp_util.launch_and_create_simple(self.repop_job_label, job_name, self.commands, command_list)

    def tutorial_contigs(self):
        print(dt.today(), "working on:", self.tutorial_string)
        assemble_contigs_path = os.path.join(self.output_folder_path, self.assemble_contigs_label)
        job_name = self.assemble_contigs_label
        command_list = self.commands.create_assemble_contigs_command(self.assemble_contigs_label, self.repop_job_label)
        self.mp_util.launch_and_create_simple(self.assemble_contigs_label, job_name, self.commands, command_list)
        mgm_file = os.path.join(assemble_contigs_path, "data", "1_mgm", "gene_report.txt")
        if(os.path.exists(mgm_file)):
            print(dt.today(), "MetaGeneMark ran successfully.  Everything is great")
        else:
            sys.exit("mgm did not run.  look into it.  Proabably a license issue. pipeline stopping here")


    def tutorial_GA(self):
        GA_BWA_start = time.time()
        GA_BWA_path = os.path.join(self.output_folder_path, self.GA_BWA_label)
        GA_BWA_jobs_folder = os.path.join(GA_BWA_path, "data", "jobs")
        #if not check_where_resume(GA_BWA_path, None, assemble_contigs_path):
        marker_path_list = []
        marker_file = "GA_split_fasta_contigs"
        marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping", marker_file)
        else:
            if(self.contig_path != "None"):
                job_name = "GA_prep_split_contigs"
                marker_path_list.append(marker_path)
                command_list = self.commands.create_split_ga_fasta_data_command(self.GA_BWA_label, self.assemble_contigs_label, "contigs", marker_file)
                self.mp_util.launch_and_create_with_mp_store(self.GA_BWA_label, job_name, self.commands, command_list)
            else:
                print(dt.today(), "no contigs supplied to GA. skipping: GA prep split contigs")
        
        sections = ["singletons"]
        if(self.read_mode == "paired"):
            sections.extend(["pair_1", "pair_2"])
        for section in sections: 
            extension = ""
            if(section == "singletons"):
                extension = os.path.splitext(self.single_path)[1]
            elif(section == "pair_1"):
                extension = os.path.splitext(self.pair_1_path)[1]
            elif(section == "pair_2"):
                extension = os.path.splitext(self.pair_2_path)[1]
            
            extension_mode = ""
            if((extension == ".fa") or (extension == ".fasta")):
                marker_file = "GA_split_fasta_" + section
                extension_mode = "fasta"
            if((extension == ".fastq") or (extension == ".fq")):
                marker_file = "GA_split_fastq_" + section
                extension_mode = "fastq"
                
            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping", marker_file)
            else:
                marker_path_list.append(marker_path)
                job_name = "GA_prep_split_" + section
                command_list = ""
                if(extension_mode == "fasta"):
                    command_list = self.commands.create_split_ga_fasta_data_command(self.GA_BWA_label, self.assemble_contigs_label, section, marker_file)
                elif(extension_mode == "fastq"):
                    command_list = self.commands.create_split_ga_fastq_data_command(self.GA_BWA_label, self.assemble_contigs_label, section, marker_file)
                self.mp_util.launch_and_create_with_mp_store(self.GA_BWA_label, job_name, self.commands, command_list)
        self.mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_prep.txt")

            
        #-------------------------------------------------------------------------
        
        sections = ["singletons"]
        if(self.contig_path != "None"):
            sections.append("contigs")
            print(dt.today(), "adding contigs to list of things to GA: BWA", sections)
            
        if self.read_mode == "paired":
            sections.extend(["pair_1", "pair_2"])
        
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                job_submitted = False
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                print("split sample:", full_sample_path)
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA" + "_" + file_tag
                marker_file = file_tag + "_bwa"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                #this checker assumes that BWA only exports a file when it's finished running
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_BWA_annotate_command_v2(self.GA_BWA_label, full_sample_path, marker_file)
                    self.mp_util.launch_and_create_with_hold(self.BWA_mem_threshold, self.BWA_job_limit, self.BWA_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)

        print(dt.today(), "all BWA jobs have launched.  waiting for them to finish")            
        self.mp_util.wait_for_mp_store()
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA.txt")
        self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            
        marker_path_list = []
        #sections = ["contigs", "singletons"]
        #if self.read_mode == "paired":
        #    sections.extend(["pair_1", "pair_2"])
        for section in sections:
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "data", "0_read_split", section)):
                full_sample_path = os.path.join(os.path.join(GA_BWA_path, "data", "0_read_split", section, split_sample))
                file_tag = os.path.basename(split_sample)
                file_tag = os.path.splitext(file_tag)[0]
                job_name = "BWA_pp" + "_" + file_tag
                marker_file = file_tag + "_bwa_pp"
                marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                    continue
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_BWA_pp_command_v2(self.GA_BWA_label, self.assemble_contigs_label, full_sample_path, marker_file)
                    self.mp_util.launch_and_create_with_hold(self.BWA_pp_mem_threshold, self.BWA_pp_job_limit, self.BWA_pp_job_delay, self.GA_BWA_label, job_name, self.commands, command_list)
                        
        print(dt.today(), "all BWA PP jobs submitted.  waiting for sync")            
        self.mp_util.wait_for_mp_store()
        if(self.contig_path != "None"):
            marker_file = "BWA_copy_contig_map"
            marker_path = os.path.join(GA_BWA_jobs_folder, marker_file)
            if(os.path.exists(marker_path)):
                print(dt.today(), "skipping:", marker_file)
            else:   
                marker_path_list.append(marker_path)
                command_list = self.commands.create_BWA_copy_contig_map_command(self.GA_BWA_label, self.assemble_contigs_label, marker_file)
                self.mp_util.launch_and_create_simple(self.GA_BWA_label, self.GA_BWA_label + "_copy_contig_map", self.commands, command_list)
        else:
            print(dt.today(), "skipping copy contig map.  There's no contigs, afterall")
        
        final_checklist = os.path.join(GA_BWA_path, "GA_BWA_pp.txt")
        self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
        
        cleanup_GA_BWA_start = time.time()
        self.mp_util.clean_or_compress(GA_BWA_path, self.keep_all, self.keep_GA_BWA)
        cleanup_GA_BWA_end = time.time()
        GA_BWA_end = time.time()
        print("GA BWA:", '%1.1f' % (GA_BWA_end - GA_BWA_start - (cleanup_GA_BWA_end - cleanup_GA_BWA_start)), "s")
        print("GA BWA cleanup:", '%1.1f' % (cleanup_GA_BWA_end - cleanup_GA_BWA_start), "s")
        
        # ------------------------------------------------
        # BLAT gene annotation
        GA_BLAT_start = time.time()
        GA_BLAT_path = os.path.join(self.output_folder_path, self.GA_BLAT_label)
        GA_BLAT_jobs_folder = os.path.join(GA_BLAT_path, "data", "jobs")
        GA_BLAT_final_job_marker = os.path.join(GA_BLAT_path, "all_BLAT")
        if (os.path.exists(GA_BLAT_final_job_marker)):
            print(dt.today(), "BLAT was run, skipping")
        else:
            
        
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    for fasta_db in os.listdir(self.paths.DNA_DB_Split):
                        if fasta_db.endswith(".fasta") or fasta_db.endswith(".ffn") or fasta_db.endswith(".fsa") or fasta_db.endswith(".fas") or fasta_db.endswith(".fna"):
                            job_name = "BLAT_" + file_tag + "_" + fasta_db
                            marker_file = file_tag + "_blat_" + fasta_db
                            marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                            #This checker assume BLAT only exports a file when it's finished running
                            if(os.path.exists(marker_path)):
                                print(dt.today(), "BLAT job ran already, skipping:", marker_file)
                                continue
                            else:
                                marker_path_list.append(marker_path)
                                command_list = self.commands.create_BLAT_annotate_command_v2(self.GA_BLAT_label, full_sample_path, fasta_db, marker_file)
                                self.mp_util.launch_only_with_hold(self.BLAT_mem_threshold, self.BLAT_job_limit, self.BLAT_job_delay, job_name, self.commands, command_list)
                                
                                    
            print(dt.today(), "final BLAT job removal")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            

                
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    job_name = file_tag + "_cat"
                    
                    marker_file = file_tag + "_blat_cat"
                    marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_BLAT_cat_command_v2(self.GA_BLAT_label, full_sample_path, marker_file)
                        self.mp_util.launch_only_with_hold(self.BLAT_mem_threshold, self.BLAT_job_limit, self.BLAT_job_delay, job_name, self.commands, command_list)
                    
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_cat.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
                
            
            
            marker_path_list = []
            for split_sample in os.listdir(os.path.join(GA_BWA_path, "final_results")):
                if(split_sample.endswith(".fasta")):
                    file_tag = os.path.basename(split_sample)
                    file_tag = os.path.splitext(file_tag)[0]
                    job_name = "BLAT_" + file_tag + "_pp"
                    full_sample_path = os.path.join(os.path.join(GA_BWA_path, "final_results", split_sample))
                    marker_file = file_tag + "_blat_pp"
                    marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                    if(os.path.exists(marker_path)):
                        print(dt.today(), "skipping:", marker_file)
                        continue
                    else:
                        marker_path_list.append(marker_path)
                        command_list = self.commands.create_BLAT_pp_command_v2(self.GA_BLAT_label, full_sample_path, self.GA_BWA_label, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.BLAT_pp_mem_threshold, self.BLAT_pp_job_limit, self.BLAT_pp_job_delay, self.GA_BLAT_label, job_name, self.commands, command_list)
                    
            print(dt.today(), "submitted all BLAT pp jobs.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            if(self.contig_path != "None"):
                job_name = "GA_BLAT_copy_contigs"
                marker_file = "blat_copy_contig_map"
                marker_path = os.path.join(GA_BLAT_jobs_folder, marker_file)
                if(os.path.exists(marker_path)):
                    print(dt.today(), "skipping:", marker_file)
                else:
                    marker_path_list.append(marker_path)
                    command_list = self.commands.create_BLAT_copy_contig_map_command(self.GA_BLAT_label, self.GA_BWA_label, marker_file)
                    self.mp_util.launch_and_create_simple(self.GA_BLAT_label, job_name, self.commands, command_list)
                
                
            final_checklist = os.path.join(GA_BLAT_path, "GA_BLAT_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            self.mp_util.delete_folder_simple(GA_BLAT_jobs_folder)
            open(GA_BLAT_final_job_marker, "a").close()
            
        cleanup_GA_BLAT_start = time.time()
        self.mp_util.clean_or_compress(GA_BLAT_path, self.keep_all, self.keep_GA_BLAT)
        cleanup_GA_BLAT_end = time.time()
        GA_BLAT_end = time.time()
        print("GA BLAT:", '%1.1f' % (GA_BLAT_end - GA_BLAT_start - (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start)), "s")
        print("GA BLAT cleanup:", '%1.1f' % (cleanup_GA_BLAT_end - cleanup_GA_BLAT_start), "s")
        
        # ------------------------------------------------------
        # Diamond gene annotation
        GA_DIAMOND_start = time.time()
        GA_DIAMOND_path = os.path.join(self.output_folder_path, self.GA_DIAMOND_label)
        GA_DIAMOND_tool_output_path = os.path.join(GA_DIAMOND_path, "data", "0_diamond")
        GA_DIAMOND_jobs_folder = os.path.join(GA_DIAMOND_path, "data", "jobs")
        
        GA_DIAMOND_final_job_marker = os.path.join(GA_DIAMOND_path, "all_dmd")
        if (os.path.exists(GA_DIAMOND_final_job_marker)):
            print(dt.today(), "DMD was run, skipping")
        else:
            
            #if not check_where_resume(None, GA_DIAMOND_tool_output_path, GA_BLAT_path, file_check_bypass = True):
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
                        command_list = self.commands.create_DIAMOND_annotate_command_v2(self.GA_DIAMOND_label, full_sample_path, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.DIAMOND_mem_threshold, self.DIAMOND_job_limit, self.DIAMOND_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                    
            print(dt.today(), "All DIAMOND jobs launched.  waiting for join")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
            open(GA_DIAMOND_final_job_marker, "a").close()
            
            
            #if not check_where_resume(GA_DIAMOND_path, None, GA_DIAMOND_tool_output_path, file_check_bypass = True):
            print(dt.today(), "DIAMOND PP threads used:", self.real_thread_count/2)
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
                        command_list = self.commands.create_DIAMOND_pp_command_v2(self.GA_DIAMOND_label, self.GA_BLAT_label, full_sample_path, marker_file)
                        self.mp_util.launch_and_create_with_hold(self.DIAMOND_pp_mem_threshold, self.DIAMOND_pp_job_limit, self.DIAMOND_pp_job_delay, self.GA_DIAMOND_label, job_name, self.commands, command_list)
                                        
            print(dt.today(), "DIAMOND pp jobs submitted.  waiting for sync")
            self.mp_util.wait_for_mp_store()
            final_checklist = os.path.join(GA_DIAMOND_path, "GA_DIAMOND_pp.txt")
            self.mp_util.check_all_job_markers(marker_path_list, final_checklist)
       
                
            
            cleanup_GA_DIAMOND_start = time.time()
            self.mp_util.clean_or_compress(GA_DIAMOND_path, self.keep_all, self.keep_GA_DIAMOND)
            cleanup_GA_DIAMOND_end = time.time()
            GA_DIAMOND_end = time.time()
            print("GA DIAMOND:", '%1.1f' % (GA_DIAMOND_end - GA_DIAMOND_start - (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start)), "s")
            print("GA DIAMOND cleanup:", '%1.1f' % (cleanup_GA_DIAMOND_end - cleanup_GA_DIAMOND_start), "s")
        
        
        GA_final_merge_start = time.time()
        GA_FINAL_MERGE_path = os.path.join(self.output_folder_path, self.GA_final_merge_label)
        marker_file = "GA_final_merge"
        marker_path = os.path.join(GA_FINAL_MERGE_path, "data", "jobs", "GA_final_merge")
        if(os.path.exists(marker_path)):
            print(dt.today(), "skipping: GA final merge")
        else:
            command_list = self.commands.create_GA_final_merge_command(self.GA_final_merge_label, self.GA_BWA_label, self.GA_BLAT_label, self.GA_DIAMOND_label, self.assemble_contigs_label, marker_file)
            job_name = "GA_final_merge"
            self.mp_util.launch_and_create_simple(self.GA_final_merge_label, job_name, self.commands, command_list)
        
        #check if all_proteins.faa was generated
        all_proteins_path = os.path.join(self.output_folder_path, self.GA_final_merge_label, "final_results", "all_proteins.faa")
        if(os.path.getsize(all_proteins_path) > 0):
            print(dt.today(), "All_proteins.faa is OK.  Continuing")
        else:
            sys.exit("GA final merge failed.  proteins weren't translated")
            
        GA_final_merge_end = time.time()
        print("GA final merge:", '%1.1f' % (GA_final_merge_end - GA_final_merge_start), "s")
        self.mp_util.clean_or_compress(GA_FINAL_MERGE_path, self.keep_all, self.keep_GA_final)