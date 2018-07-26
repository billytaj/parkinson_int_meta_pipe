#The functions here generate the pipeline commands.
#Each command module is made up of (often) many sub stages that are used to get the final result.
#If you want to move around the ordering, you'd do that here.

import os
import sys
import mt_pipe_paths as mpp
import subprocess as sp
#------------------------------------------------------
# pragmas needed for command construction

SLURM_Submit = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=4:00:00
#SBATCH --job-name=NAME
#SBATCH --error=ERROR
#SBATCH --output=OUTPUT
 
cd $SLURM_SUBMIT_DIR
 
module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras anaconda3/4.0.0
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

OLDPATH=$PATH:/home/j/jparkin/ctorma/emboss/bin/:/home/j/jparkin/mobolaji/Tools/Barrnap/bin/:/home/j/jparkin/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkin/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkin/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkin/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

#--------------------------------------------------------------------------
# class definition

class mt_pipe_commands:
    #--------------------------------------------------------------------
    # constructor:
    # there should only be one of these objects used for an entire pipeline.
    def __init__(self, Quality_score = 33, Thread_count = 80, system_mode = "scinet", user_mode = "billy",  raw_sequence_path_0 = None, raw_sequence_path_1 = None):
        # path to the raw genome sequence file
        Input_File = os.getcwd()
        if not(raw_sequence_path_0 is None):
            self.raw_sequence_path_0 = raw_sequence_path_0
            print("raw sequence 0:", self.raw_sequence_path_0)

        if not(raw_sequence_path_1 is None):
            self.raw_sequence_path_1 = raw_sequence_path_1
            print("raw seqeunce 1:", self.raw_sequence_path_1)
        self.system_mode = system_mode
        self.tool_path_obj = mpp.tool_path_obj(system_mode, user_mode)
        
        self.Qual_str = str(Quality_score)
        self.Input_Path = os.path.dirname(Input_File)
        self.Threads_str= str(Thread_count)
        
        print("input filepath:", Input_File)
    
    #-----------------------------------------------------------
    # support functions
    #def make_stage_dir(self, stage_name):
    #    if not(os.path.exists(self.)):
    def make_folder(self, folder_path):
        if not(os.path.exists(folder_path)):
            os.makedirs(folder_path)

    def create_pbs_and_launch(self, job_name, command_list, run_job = False, inner_name = None, mode = "low", dependency_list = None,  work_in_background = False):
        #create the pbs job, and launch items
        #job name: string tag for export file name
        #command list:  list of command statements for writing
        #mode: selection of which pbs template to use: default -> low memory
        #dependency_list: if not empty, will append wait args to sbatch subprocess call. it's polymorphic
        #returns back the job ID given from sbatch

        if self.system_mode == "scinet":
            pbs_template = SLURM_Submit

        #    if(mode == "med"):
        #        pbs_template = PBS_Submit_HighMem
        #    elif(mode == "high"):
        #        pbs_template = PBS_Submit_vHighMem
        #    else:
        #        pbs_template = PBS_Submit_LowMem

        #job name:              string tag for export file name
        #command list:          list of command statements for writing
        #run_job:               the ability to just generate the shell, and not run it
        #inner_name:            to make it so that a new shellscript is generated for each split job -> name override
        #work_in_background:    the ability to run the job in background.  was used in single-cpu mode, but no longer needed
        
        #returns nothing

        #if(self.system_mode == "singularity" or self.system_mode == "docker"):
            #docker mode, depending on the machine it's applied to, it may be multi-core
            pbs_script_full_path = os.getcwd() + "/" + job_name +"/" + job_name
            if not inner_name is None:
                pbs_script_full_path = os.getcwd() + "/" + job_name + "/" + inner_name

            try:
                with open(pbs_script_full_path + ".sh", "w+") as PBS_script_out:
                    for line in pbs_template.splitlines():
                        if "NAME" in line:
                            line = line.replace("NAME", pbs_script_full_path)
                        if "ERROR" in line:
                            line = line.replace("ERROR", pbs_script_full_path + "_ERR")
                        if "OUTPUT" in line:
                            line = line.replace("OUTPUT", pbs_script_full_path + "_OUT")
                        if "COMMANDS" in line:
                            PBS_script_out.write("\n".join(command_list))
                            break

                        PBS_script_out.write(line + "\n")
                    PBS_script_out.close()
                    dep_str = ""
                    if isinstance(dependency_list, int):
                        #single dep
                        dep_str = "-W depend=afterok:" + str(dependency_list)
                        print(dep_str)
                        if inner_name is None:
                            print(job_name, "running with single dependency")
                        else:
                            print(inner_name, "running with single dependency")

                    elif isinstance(dependency_list, list):
                        # multiple deps
                        print("mutiple dependencies being used")
                        dep_str = "-d depend=afterok"
                        for item in dependency_list:
                            dep_str += ":"+str(item)
                        if inner_name is None:
                            print(job_name, "running with multiple dependencies")
                        else:
                            print(inner_name, "running with multiple dependencies")
                    elif dependency_list is None:
                        if inner_name is None:
                            print(job_name, "running without dependency")
                        else:
                            print(inner_name, "running without dependency")
                    else:
                        print("This isn't supposed to happen")

                    if run_job: # a lock built for testing syntax, but not run
                        try:
                            if not dep_str == "":
                                print("dep string not empty")
                                job_id = sp.check_output(["sbatch", pbs_script_full_path+".sh", dep_str])
                            else:
                                job_id = sp.check_output(["sbatch", pbs_script_full_path+".sh"])
                            #return val is a binary string, so we convert, and extract only the numeric part
                            job_id = job_id.decode('ascii')
                            job_id = int(job_id.split(' ')[-1])

                            return job_id
                        except Exception as e:
                            print("subprocess call error:", e)
                    else:
                        return 0


            except Exception as e:
                # error catchall
                print("Failure at pbs creation:", e)
                sys.exit()
        else:
            #docker mode: single cpu
            # no ID, no sbatch.  just run the command
            pbs_script_full_path = os.getcwd() + "/" + job_name +"/" + job_name
            if not inner_name is None:
                pbs_script_full_path = os.getcwd() + "/" + job_name + "/" + inner_name
            with open(pbs_script_full_path + ".sh", "w+") as PBS_script_out:
                for item in command_list:
                    PBS_script_out.write(item + "\n")
                PBS_script_out.close()
            if run_job:
                if not work_in_background:
                    try:
                        sp.check_output(["sh", pbs_script_full_path + ".sh"])
                    except sp.CalledProcessError as e:
                        return_code = e.returncode
                        if return_code != 1:
                            raise
                else:
                    try:
                        process_id = sp.Popen(["sh", pbs_script_full_path + ".sh"])
                        return process_id
                    except sp.CalledProcessError as e:
                        return_code = e.returncode
                        if return_code != 1:
                            raise
            else:
                print("not running job.  run_job set to False")

    def create_pre_double_command(self, stage_name):
        #why do we leave all the interim files intact?
        #because science needs repeatable data, and the process needs to be able to start at any point
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        if not (os.path.exists(subfolder)):
            os.makedirs(subfolder)
        if not(os.path.exists(data_folder)):
            os.makedirs(data_folder)

        sorted_raw_folder = data_folder + "0_sorted_raw_input/"
        self.make_folder(sorted_raw_folder)

        sort_pair_1 = ">&2 echo Sorting pair 1 | "
        sort_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.sort_reads + " "
        sort_pair_1 += self.raw_sequence_path_0
        sort_pair_1 += " " + sorted_raw_folder + "pair_1_sorted.fastq"

        sort_pair_2 = ">&2 echo Sorting pair 2 | "
        sort_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.sort_reads + " "
        sort_pair_2 += self.raw_sequence_path_1
        sort_pair_2 += " " + sorted_raw_folder + "pair_2_sorted.fastq"

        adapter_folder = data_folder + "1_adapter_removal/"
        self.make_folder(adapter_folder)
        adapter_removal_line = ">&2 echo Removing adapters | "
        adapter_removal_line += self.tool_path_obj.AdapterRemoval
        #adapter_removal_line += " --file1 " + self.raw_sequence_path_0
        #adapter_removal_line += " --file2 " + self.raw_sequence_path_1
        adapter_removal_line += " --file1 " + sorted_raw_folder + "pair_1_sorted.fastq"
        adapter_removal_line += " --file2 " + sorted_raw_folder + "pair_2_sorted.fastq"
        adapter_removal_line += " --qualitybase " + str(self.Qual_str) #must be either 33 or 64
        adapter_removal_line += " --threads " + self.Threads_str
        adapter_removal_line += " --minlength " + "30"
        adapter_removal_line += " --basename " + adapter_folder #+ os.path.splitext(self.Input_FName)[0]
        adapter_removal_line += "_AdapterRemoval"
        adapter_removal_line += " --trimqualities "
        adapter_removal_line += " --output1 " + adapter_folder + "pair_1_adptr_rem.fastq"
        adapter_removal_line += " --output2 " + adapter_folder + "pair_2_adptr_rem.fastq"
        adapter_removal_line += " --singleton " + adapter_folder + "singletons_adptr_rem.fastq"

        # tries to merge the cleaned pairs
        # rejects get sent out
        vsearch_merge_folder = data_folder + "2_vsearch_pair_merge/"
        self.make_folder(vsearch_merge_folder)
        vsearch_merge = ">&2 echo " + "Vsearch Merge pairs | "
        vsearch_merge += self.tool_path_obj.vsearch
        vsearch_merge += " --fastq_mergepairs " + adapter_folder + "pair_1_adptr_rem.fastq"
        vsearch_merge += " --reverse " + adapter_folder + "pair_2_adptr_rem.fastq"
        vsearch_merge += " --fastq_ascii " + str(self.Qual_str)
        vsearch_merge += " --fastqout " + vsearch_merge_folder + "merge_success.fastq"
        vsearch_merge += " --fastqout_notmerged_fwd " + vsearch_merge_folder + "pair_1_merge_reject.fastq"
        vsearch_merge += " --fastqout_notmerged_rev " + vsearch_merge_folder + "pair_2_merge_reject.fastq"

        # concatenate the merge overlaps with the singletons
        cat_glue = ">&2 echo concatenating orphans | "
        cat_glue += "cat "
        cat_glue += vsearch_merge_folder + "merge_success.fastq "
        cat_glue += adapter_folder + "singletons_adptr_rem.fastq"
        cat_glue += " > " + vsearch_merge_folder + "orphans.fastq"


        #Filter out low-quality reads
        #start with the singles / merged sections
        vsearch_filter_folder = data_folder + "3_ar_quality_filter/"
        self.make_folder(vsearch_filter_folder)
        vsearch_filter_0 = ">&2 echo low-quality filter on orphans | "
        vsearch_filter_0 += self.tool_path_obj.vsearch
        vsearch_filter_0 += " --fastq_filter " + vsearch_merge_folder + "orphans.fastq"
        vsearch_filter_0 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_0 += " --fastq_maxee " + "2.0"
        vsearch_filter_0 += " --fastqout " + vsearch_filter_folder + "orphans_hq.fastq"

        #then move onto the standalones in pair 1
        vsearch_filter_1 = ">&2 echo low-quality filter on pair 1 | "
        vsearch_filter_1 += self.tool_path_obj.vsearch
        vsearch_filter_1 += " --fastq_filter " + vsearch_merge_folder + "pair_1_merge_reject.fastq"
        vsearch_filter_1 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_1 += " --fastq_maxee " + "2.0"
        vsearch_filter_1 += " --fastqout " + vsearch_filter_folder + "pair_1_hq.fastq"

        vsearch_filter_2 = ">&2 echo low-quality filter on pair 2 | "
        vsearch_filter_2 += self.tool_path_obj.vsearch
        vsearch_filter_2 += " --fastq_filter " + vsearch_merge_folder + "pair_2_merge_reject.fastq"
        vsearch_filter_2 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_2 += " --fastq_maxee " + "2.0"
        vsearch_filter_2 += " --fastqout " + vsearch_filter_folder + "pair_2_hq.fastq"

        #redistribute data into orphans, or paired-reads
        orphan_read_filter_folder = data_folder + "4_orphan_read_filter/"
        self.make_folder(orphan_read_filter_folder)
        orphan_read_filter = ">&2 echo moving newly orphaned hq reads | "
        orphan_read_filter += self.tool_path_obj.Python + " "
        orphan_read_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_read_filter += vsearch_filter_folder + "pair_1_hq.fastq "
        orphan_read_filter += vsearch_filter_folder + "pair_2_hq.fastq "
        orphan_read_filter += vsearch_filter_folder + "orphans_hq.fastq "
        orphan_read_filter += orphan_read_filter_folder + "pair_1_match.fastq "
        orphan_read_filter += orphan_read_filter_folder + "pair_2_match.fastq "
        orphan_read_filter += orphan_read_filter_folder + "orphans.fastq"

        #remove duplicates (to shrink the data size)
        cdhit_folder = data_folder + "5_remove_duplicates/"
        self.make_folder(cdhit_folder)
        cdhit_orphans = ">&2 echo removing orphan duplicates | "
        cdhit_orphans += self.tool_path_obj.cdhit_dup + " -i "
        cdhit_orphans += orphan_read_filter_folder + "orphans.fastq"
        cdhit_orphans += " -o " + cdhit_folder + "orphans_unique.fastq"


        #move_unpaired_cluster = "mv "
        #move_unpaired_cluster += self.Input_Filepath + "_unpaired_unique.fastq.clstr "
        #move_unpaired_cluster += self.Input_Filepath + "_unpaired.clstr"

        #remove duplicates in the pairs
        cdhit_pair_1 = ">&2 echo remove duplicates from pair 1 | "
        cdhit_pair_1 += self.tool_path_obj.cdhit_dup
        cdhit_pair_1 += " -i " + orphan_read_filter_folder + "pair_1_match.fastq"
        cdhit_pair_1 += " -o " + cdhit_folder + "pair_1_unique.fastq"

        cdhit_pair_2 = ">&2 echo remove duplicates from pair 2 | "
        cdhit_pair_2 += self.tool_path_obj.cdhit_dup
        cdhit_pair_2 += " -i " + orphan_read_filter_folder + "pair_2_match.fastq"
        cdhit_pair_2 += " -o " + cdhit_folder + "pair_2_unique.fastq"

        #move_paired_cluster = "mv " + self.Input_File1 + "_paired_unique.fastq.clstr" + " " + self.Input_Filepath + "_paired.clstr"
        #NOTE: we overwrite host contaminants here. why? because the preprocess folder is created
        #in this function only
        host_removal_folder = data_folder + "6_remove_host/"
        self.make_folder(host_removal_folder)

        Host_Contaminants = host_removal_folder + "host_contaminents_seq.fasta"
        copy_host = ">&2 echo Copy the host file over | "
        copy_host += "cp " + self.tool_path_obj.Host + " " + Host_Contaminants

        #craft a BWA index for the host sequences
        bwa_hr_prep = ">&2 echo make host contaminants index for BWA | "
        bwa_hr_prep += self.tool_path_obj.BWA + " index -a bwtsw " + Host_Contaminants

        samtools_hr_prep = ">&2 echo SAMTOOLS host contaminant prep | "
        samtools_hr_prep += self.tool_path_obj.SAMTOOLS + " faidx " + Host_Contaminants

        #host removal on unique orphans
        bwa_hr_orphans = ">&2 echo BWA host remove on orphans | "
        bwa_hr_orphans += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_orphans += self.Threads_str + " "
        bwa_hr_orphans += Host_Contaminants + " "
        bwa_hr_orphans += cdhit_folder + "orphans_unique.fastq"
        bwa_hr_orphans += " > " + host_removal_folder + "orphans_no_host.sam"

        #annoying type conversion pt 1
        samtools_hr_orphans_sam_to_bam = ">&2 echo convert orphans hr files pt1 | "
        samtools_hr_orphans_sam_to_bam += self.tool_path_obj.SAMTOOLS
        samtools_hr_orphans_sam_to_bam += " view -bS " + host_removal_folder + "orphans_no_host.sam"
        samtools_hr_orphans_sam_to_bam += " > " + host_removal_folder + "orphans_no_host.bam"
        #annoying type conversion pt 2
        samtools_no_host_orphans_bam_to_fastq = ">&2 echo convert orphans hr files pt2 | "
        samtools_no_host_orphans_bam_to_fastq += self.tool_path_obj.SAMTOOLS
        samtools_no_host_orphans_bam_to_fastq += " fastq -n -f 4" + " -0 " + host_removal_folder + "orphans_no_host.fastq" + " "
        samtools_no_host_orphans_bam_to_fastq += host_removal_folder + "orphans_no_host.bam"

        #apparently, we're to keep the host separation
        samtools_host_orphans_bam_to_fastq = ">&2 echo convert orphans hr files pt3 | "
        samtools_host_orphans_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_host_orphans_bam_to_fastq += " -0 " + host_removal_folder + "orphans_host_only.fastq" + " "
        samtools_host_orphans_bam_to_fastq += host_removal_folder + "orphans_no_host.bam"

        #host-remove the rest
        bwa_hr_pair = ">&2 echo bwa pair host remove | "
        bwa_hr_pair += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_pair += self.Threads_str + " "
        bwa_hr_pair += Host_Contaminants + " "
        bwa_hr_pair += cdhit_folder + "pair_1_unique.fastq" + " "
        bwa_hr_pair += cdhit_folder + "pair_2_unique.fastq"
        bwa_hr_pair += " > " + host_removal_folder + "pair_no_host.sam"

        # bwa hr pair 1 only
        bwa_hr_pair_1 = ">&2 echo bwa pair host remove | "
        bwa_hr_pair_1 += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_pair_1 += self.Threads_str + " "
        bwa_hr_pair_1 += Host_Contaminants + " "
        bwa_hr_pair_1 += cdhit_folder + "pair_1_unique.fastq"
        bwa_hr_pair_1 += " > " + host_removal_folder + "pair_1_no_host.sam"


        #separating bwa results back into paired reads
        samtools_host_pair_1_sam_to_bam = ">&2 echo convert pair hr files pt1 | "
        samtools_host_pair_1_sam_to_bam += self.tool_path_obj.SAMTOOLS + " view -bS " + host_removal_folder + "pair_1_no_host.sam"
        samtools_host_pair_1_sam_to_bam += " > " + host_removal_folder + "pair_1_no_host.bam"

        #stuff that doesn't match with the host
        samtools_no_host_pair_1_bam_to_fastq = ">&2 echo convert pair hr files pt2 | "
        samtools_no_host_pair_1_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_host_pair_1_bam_to_fastq += " -0 " + host_removal_folder + "pair_1_no_host.fastq" # out
        samtools_no_host_pair_1_bam_to_fastq += " " + host_removal_folder + "pair_1_no_host.bam" #in

        #stuff that matches with the host (why keep it?  request from john)
        samtools_host_pair_1_bam_to_fastq = ">&2 echo convert pair hr files pt3 | "
        samtools_host_pair_1_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_host_pair_1_bam_to_fastq += " -0 " + host_removal_folder + "pair_1_host_only.fastq"
        samtools_host_pair_1_bam_to_fastq += " " + host_removal_folder + "pair_1_no_host.bam"


        # bwa hr pair 1 only
        bwa_hr_pair_2 = ">&2 echo bwa pair host remove | "
        bwa_hr_pair_2 += self.tool_path_obj.BWA + " mem -t "
        bwa_hr_pair_2 += self.Threads_str + " "
        bwa_hr_pair_2 += Host_Contaminants + " "
        bwa_hr_pair_2 += cdhit_folder + "pair_2_unique.fastq"
        bwa_hr_pair_2 += " > " + host_removal_folder + "pair_2_no_host.sam"

        #separating bwa results back into paired reads
        samtools_host_pair_2_sam_to_bam = ">&2 echo convert pair hr files pt1 | "
        samtools_host_pair_2_sam_to_bam += self.tool_path_obj.SAMTOOLS + " view -bS " + host_removal_folder + "pair_2_no_host.sam"
        samtools_host_pair_2_sam_to_bam += " > " + host_removal_folder + "pair_2_no_host.bam"

        #stuff that doesn't match with the host
        samtools_no_host_pair_2_bam_to_fastq = ">&2 echo convert pair hr files pt2 | "
        samtools_no_host_pair_2_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_host_pair_2_bam_to_fastq += " -0 " + host_removal_folder + "pair_2_no_host.fastq" # out
        samtools_no_host_pair_2_bam_to_fastq += " " + host_removal_folder + "pair_2_no_host.bam" #in

        #stuff that matches with the host (why keep it?  request from john)
        samtools_host_pair_2_bam_to_fastq = ">&2 echo convert pair hr files pt3 | "
        samtools_host_pair_2_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_host_pair_2_bam_to_fastq += " -0 " + host_removal_folder + "pair_2_host_only.fastq"
        samtools_host_pair_2_bam_to_fastq += " " + host_removal_folder + "pair_2_no_host.bam"

        #blat prep
        make_blast_db_host = ">&2 echo Make BLAST db for host contaminants | "
        make_blast_db_host += self.tool_path_obj.Makeblastdb + " -in " + Host_Contaminants + " -dbtype nucl"

        vsearch_filter_3 = ">&2 echo Convert orphans for BLAT | "
        vsearch_filter_3 += self.tool_path_obj.vsearch
        vsearch_filter_3 += " --fastq_filter " + host_removal_folder + "orphans_no_host.fastq"
        vsearch_filter_3 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_3 += " --fastaout " + host_removal_folder + "orphans_no_host.fasta"

        vsearch_filter_4 = ">&2 echo Convert pair 1 for BLAT | "
        vsearch_filter_4 += self.tool_path_obj.vsearch
        vsearch_filter_4 += " --fastq_filter " + host_removal_folder + "pair_1_no_host.fastq"
        vsearch_filter_4 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_4 += " --fastaout " + host_removal_folder + "pair_1_no_host.fasta"

        vsearch_filter_5 = ">&2 echo Convert pair 2 for BLAT | "
        vsearch_filter_5 += self.tool_path_obj.vsearch
        vsearch_filter_5 += " --fastq_filter " + host_removal_folder + "pair_2_no_host.fastq"
        vsearch_filter_5 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_5 += " --fastaout " + host_removal_folder + "pair_2_no_host.fasta"

        blat_hr_orphans = ">&2 echo BLAT hr orphans | "
        blat_hr_orphans += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_hr_orphans += Host_Contaminants + " "
        blat_hr_orphans += host_removal_folder + "orphans_no_host.fasta"
        blat_hr_orphans += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_orphans += " " + host_removal_folder + "orphans_no_host.blatout"

        blat_hr_pair_1 = ">&2 echo BLAT hr pair 1 | "
        blat_hr_pair_1 += self.tool_path_obj.BLAT
        blat_hr_pair_1 += " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " "
        blat_hr_pair_1 += host_removal_folder + "pair_1_no_host.fasta"
        blat_hr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_1 += " " + host_removal_folder + "pair_1_no_host.blatout"

        blat_hr_pair_2 = ">&2 echo BLAT hr pair 2 | "
        blat_hr_pair_2 += self.tool_path_obj.BLAT
        blat_hr_pair_2 += " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " "
        blat_hr_pair_2 += host_removal_folder + "pair_2_no_host.fasta"
        blat_hr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_2 += " " + host_removal_folder + "pair_2_no_host.blatout"

        # HR BLAT
        blat_hr_folder = data_folder + "7_blat_hr/"
        self.make_folder(blat_hr_folder)
        hr_orphans = ">&2 echo BLAT contaminant orphans | "
        hr_orphans += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        hr_orphans += host_removal_folder + "orphans_no_host.fastq" + " " # in
        hr_orphans += host_removal_folder + "orphans_no_host.blatout" + " " #in
        hr_orphans += blat_hr_folder + "orphans_no_host.fastq" + " " #out
        hr_orphans += blat_hr_folder + "orphans_host_only.fastq" #out

        hr_pair_1 = ">&2 echo BLAT contaminant pair 1 | "
        hr_pair_1 += self.tool_path_obj.Python + " "
        hr_pair_1 += self.tool_path_obj.BLAT_Contaminant_Filter + " "
        hr_pair_1 += host_removal_folder + "pair_1_no_host.fastq" + " "
        hr_pair_1 += host_removal_folder + "pair_1_no_host.blatout" + " "
        hr_pair_1 += blat_hr_folder + "pair_1_no_host.fastq" + " "
        hr_pair_1 += blat_hr_folder + "pair_1_host_only.fastq"

        hr_pair_2 = ">&2 echo BLAT contaminant pair 2 | "
        hr_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        hr_pair_2 += host_removal_folder + "pair_2_no_host.fastq" + " "
        hr_pair_2 += host_removal_folder + "pair_2_no_host.blatout" + " "
        hr_pair_2 += blat_hr_folder + "pair_2_no_host.fastq" + " "
        hr_pair_2 += blat_hr_folder + "pair_2_host_only.fastq"

        #-----------------------------
        # Vector removal
        # vectors are artifacts that came from the read-gathering
        vector_removal_folder = data_folder + "8_vector_removal/"
        self.make_folder(vector_removal_folder)
        Vector_Contaminants = vector_removal_folder + "vector_contaminants_seq.fasta"

        copy_vector = ">&2 echo copy vector prep | "
        copy_vector += "cp " + self.tool_path_obj.UniVec_Core + " " + Vector_Contaminants

        bwa_vr_prep = ">&2 echo BWA vector prep | "
        bwa_vr_prep += self.tool_path_obj.BWA + " index -a bwtsw " + Vector_Contaminants

        samtools_vr_prep = ">&2 echo samtools vector prep | "
        samtools_vr_prep += self.tool_path_obj.SAMTOOLS + " faidx " + Vector_Contaminants

        bwa_vr_orphans = ">&2 echo BWA vr oprhans | "
        bwa_vr_orphans += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_orphans += Vector_Contaminants + " "
        bwa_vr_orphans += blat_hr_folder + "orphans_no_host.fastq"
        bwa_vr_orphans += " > " + vector_removal_folder + "orphans_no_vectors.sam"

        samtools_no_vector_orphans_sam_to_bam = ">&2 echo samtools vr oprhans pt 1 | "
        samtools_no_vector_orphans_sam_to_bam += self.tool_path_obj.SAMTOOLS + " view -bS "
        samtools_no_vector_orphans_sam_to_bam += vector_removal_folder + "orphans_no_vectors.sam"
        samtools_no_vector_orphans_sam_to_bam += " > " + vector_removal_folder + "orphans_no_vectors.bam"

        samtools_no_vector_orphans_bam_to_fastq = ">&2 echo samtools vr orphans pt 2 | "
        samtools_no_vector_orphans_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_vector_orphans_bam_to_fastq += " -0 " + vector_removal_folder + "orphans_no_vectors.fastq "
        samtools_no_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_no_vectors.bam"

        samtools_vector_orphans_bam_to_fastq = ">&2 echo samtools vr orphans pt 3 | "
        samtools_vector_orphans_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_vector_orphans_bam_to_fastq += " -0 " + vector_removal_folder + "orphans_vectors_only.fastq "
        samtools_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_no_vectors.bam"

        bwa_vr_pair_1 = ">&2 echo bwa vr pair | "
        bwa_vr_pair_1 += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_pair_1 += Vector_Contaminants + " "
        bwa_vr_pair_1 += blat_hr_folder + "pair_1_no_host.fastq "
        bwa_vr_pair_1 += " > " + vector_removal_folder + "pair_1_no_vectors.sam"

        bwa_vr_pair_2 = ">&2 echo bwa vr pair | "
        bwa_vr_pair_2 += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " "
        bwa_vr_pair_2 += Vector_Contaminants + " "
        bwa_vr_pair_2 += blat_hr_folder + "pair_2_no_host.fastq "
        bwa_vr_pair_2 += " > " + vector_removal_folder + "pair_2_no_vectors.sam"


        samtools_vr_pair_1_sam_to_bam = ">&2 echo samtools vr pair pt 1 | "
        samtools_vr_pair_1_sam_to_bam += self.tool_path_obj.SAMTOOLS + " view -bS "
        samtools_vr_pair_1_sam_to_bam += vector_removal_folder + "pair_1_no_vectors.sam"
        samtools_vr_pair_1_sam_to_bam += " > " + vector_removal_folder + "pair_1_no_vectors.bam"

        samtools_no_vector_pair_1_bam_to_fastq = ">&2 echo samtools vr pair pt 2 | "
        samtools_no_vector_pair_1_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_vector_pair_1_bam_to_fastq += " -0 " + vector_removal_folder + "pair_1_no_vectors.fastq"
        samtools_no_vector_pair_1_bam_to_fastq += " " + vector_removal_folder + "pair_1_no_vectors.bam"

        samtools_vector_pair_1_bam_to_fastq = ">&2 echo samtools vr pair pt 3 | "
        samtools_vector_pair_1_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_vector_pair_1_bam_to_fastq += " -0 " + vector_removal_folder + "pair_1_vectors_only.fastq"
        samtools_vector_pair_1_bam_to_fastq += " " + vector_removal_folder + "pair_1_no_vectors.bam"

        samtools_vr_pair_2_sam_to_bam = ">&2 echo samtools vr pair pt 1 | "
        samtools_vr_pair_2_sam_to_bam += self.tool_path_obj.SAMTOOLS + " view -bS "
        samtools_vr_pair_2_sam_to_bam += vector_removal_folder + "pair_2_no_vectors.sam"
        samtools_vr_pair_2_sam_to_bam += " > " + vector_removal_folder + "pair_2_no_vectors.bam"

        samtools_no_vector_pair_2_bam_to_fastq = ">&2 echo samtools vr pair pt 2 | "
        samtools_no_vector_pair_2_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -f 4"
        samtools_no_vector_pair_2_bam_to_fastq += " -0 " + vector_removal_folder + "pair_2_no_vectors.fastq"
        samtools_no_vector_pair_2_bam_to_fastq += " " + vector_removal_folder + "pair_2_no_vectors.bam"

        samtools_vector_pair_2_bam_to_fastq = ">&2 echo samtools vr pair pt 3 | "
        samtools_vector_pair_2_bam_to_fastq += self.tool_path_obj.SAMTOOLS + " fastq -n -F 4"
        samtools_vector_pair_2_bam_to_fastq += " -0 " + vector_removal_folder + "pair_2_vectors_only.fastq"
        samtools_vector_pair_2_bam_to_fastq += " " + vector_removal_folder + "pair_2_no_vectors.bam"



        make_blast_db_vector = ">&2 echo BLAST make db vectors | "
        make_blast_db_vector += self.tool_path_obj.Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl"

        vsearch_filter_6 = ">&2 echo convert vr orphans for BLAT | "
        vsearch_filter_6 += self.tool_path_obj.vsearch
        vsearch_filter_6 += " --fastq_filter " + vector_removal_folder + "orphans_no_vectors.fastq"
        vsearch_filter_6 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_6 += " --fastaout " + vector_removal_folder + "orphans_no_vectors.fasta"

        vsearch_filter_7 = ">&2 echo convert vr pair 1 for BLAT | "
        vsearch_filter_7 += self.tool_path_obj.vsearch
        vsearch_filter_7 += " --fastq_filter " + vector_removal_folder + "pair_1_no_vectors.fastq"
        vsearch_filter_7 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_7 += " --fastaout " + vector_removal_folder + "pair_1_no_vectors.fasta"

        vsearch_filter_8 = ">&2 echo convert vr pair 2 for BLAT | "
        vsearch_filter_8 += self.tool_path_obj.vsearch
        vsearch_filter_8 += " --fastq_filter " + vector_removal_folder + "pair_2_no_vectors.fastq"
        vsearch_filter_8 += " --fastq_ascii " + self.Qual_str
        vsearch_filter_8 += " --fastaout " + vector_removal_folder + "pair_2_no_vectors.fasta"

        blat_vr_orphans = ">&2 echo BLAT vr orphans | "
        blat_vr_orphans += self.tool_path_obj.BLAT
        blat_vr_orphans += " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_orphans += Vector_Contaminants + " "
        blat_vr_orphans += vector_removal_folder + "orphans_no_vectors.fasta"
        blat_vr_orphans += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_orphans += vector_removal_folder + "orphans_no_vectors.blatout"

        blat_vr_pair_1 = ">&2 echo BLAT vr pair 1 | "
        blat_vr_pair_1 += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_1 += Vector_Contaminants + " "
        blat_vr_pair_1 += vector_removal_folder + "pair_1_no_vectors.fasta"
        blat_vr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_pair_1 += vector_removal_folder + "pair_1_no_vectors.blatout"

        blat_vr_pair_2 = ">&2 echo BLAT vr pair 2 | "
        blat_vr_pair_2 += self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_vr_pair_2 += Vector_Contaminants + " "
        blat_vr_pair_2 += vector_removal_folder + "pair_2_no_vectors.fasta"
        blat_vr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " "
        blat_vr_pair_2 += vector_removal_folder + "pair_2_no_vectors.blatout"

        blat_containment_vector_folder = data_folder + "9_blat_containment_vr/"
        self.make_folder(blat_containment_vector_folder)
        blat_containment_vector_orphans = ">&2 echo BLAT contaminant orphans | "
        blat_containment_vector_orphans += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        blat_containment_vector_orphans += vector_removal_folder + "orphans_no_vectors.fastq" + " " #in
        blat_containment_vector_orphans += vector_removal_folder + "orphans_no_vectors.blatout" + " " #in
        blat_containment_vector_orphans += blat_containment_vector_folder + "orphans_no_vectors.fastq" + " " #out
        blat_containment_vector_orphans += blat_containment_vector_folder + "orphans_vectors_only.fastq" #out

        blat_containment_vector_pair_1 = ">&2 echo BLAT contaminant pair 1 | "
        blat_containment_vector_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        blat_containment_vector_pair_1 += vector_removal_folder + "pair_1_no_vectors.fastq" + " "
        blat_containment_vector_pair_1 += vector_removal_folder + "pair_1_no_vectors.blatout" + " "
        blat_containment_vector_pair_1 += blat_containment_vector_folder + "pair_1_no_vectors.fastq" + " "
        blat_containment_vector_pair_1 += blat_containment_vector_folder + "pair_1_vectors_only.fastq"

        blat_containment_vector_pair_2 = ">&2 echo BLAT contaminant pair 2 | "
        blat_containment_vector_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.BLAT_Contaminant_Filter + " "
        blat_containment_vector_pair_2 += vector_removal_folder + "pair_2_no_vectors.fastq" + " "
        blat_containment_vector_pair_2 += vector_removal_folder + "pair_2_no_vectors.blatout" + " "
        blat_containment_vector_pair_2 += blat_containment_vector_folder + "pair_2_no_vectors.fastq" + " "
        blat_containment_vector_pair_2 += blat_containment_vector_folder + "pair_2_vectors_only.fastq"

        #----final split
        #----------------------------------
        #our convention dictates that all final results must be placed inside this folder.
        final_folder = data_folder + "final_results/"
        self.make_folder(final_folder)

        move_orphans = "mv " + blat_containment_vector_folder + "orphans_no_vectors.fastq "
        move_orphans += final_folder

        move_pair_1 = "mv " + blat_containment_vector_folder + "pair_1_no_vectors.fastq "
        move_pair_1 += final_folder

        move_pair_2 = "mv " + blat_containment_vector_folder + "pair_2_no_vectors.fastq "
        move_pair_2 += final_folder



        COMMANDS_Pre = [
            #make the pairs align, by sorting
            sort_pair_1,
            sort_pair_2,
            # # remove adapters
            adapter_removal_line,
            # #trim things
            vsearch_merge,
            cat_glue,
            vsearch_filter_0,
            vsearch_filter_1,
            vsearch_filter_2,
            orphan_read_filter,

            cdhit_orphans,
            # # # move_unpaired_cluster,
            cdhit_pair_1,
            cdhit_pair_2,
            # # # move_paired_cluster,
            # # #----host removal
            copy_host,
            bwa_hr_prep,
            # # #----SAMTOOLS makes bam files
            samtools_hr_prep,
            bwa_hr_orphans,
            samtools_hr_orphans_sam_to_bam,
            samtools_no_host_orphans_bam_to_fastq,
            samtools_host_orphans_bam_to_fastq,
            bwa_hr_pair_1,
            samtools_host_pair_1_sam_to_bam,
            samtools_no_host_pair_1_bam_to_fastq,
            samtools_host_pair_1_bam_to_fastq,
            bwa_hr_pair_2,
            samtools_host_pair_2_sam_to_bam,
            samtools_no_host_pair_2_bam_to_fastq,
            samtools_host_pair_2_bam_to_fastq,
            make_blast_db_host,
            vsearch_filter_3,
            vsearch_filter_4,
            vsearch_filter_5,
            blat_hr_orphans,
            blat_hr_pair_1,
            blat_hr_pair_2,
            hr_orphans,
            hr_pair_1,
            hr_pair_2,
            # # # #-----vector removal
            copy_vector,
            bwa_vr_prep,
            samtools_vr_prep,
            bwa_vr_orphans,
            samtools_no_vector_orphans_sam_to_bam,
            samtools_no_vector_orphans_bam_to_fastq,
            samtools_vector_orphans_bam_to_fastq,
            bwa_vr_pair_1,
            bwa_vr_pair_2,
            samtools_vr_pair_1_sam_to_bam,
            samtools_no_vector_pair_1_bam_to_fastq,
            samtools_vector_pair_1_bam_to_fastq,
            samtools_vr_pair_2_sam_to_bam,
            samtools_no_vector_pair_2_bam_to_fastq,
            samtools_vector_pair_2_bam_to_fastq,

            make_blast_db_vector,
            vsearch_filter_6,
            vsearch_filter_7,
            vsearch_filter_8,
            blat_vr_orphans,
            blat_vr_pair_1,
            blat_vr_pair_2,
            blat_containment_vector_orphans,
            blat_containment_vector_pair_1,
            blat_containment_vector_pair_2,
            move_orphans,
            move_pair_1,
            move_pair_2

        ]

        return COMMANDS_Pre

    def create_rRNA_filter_prep_command(self, stage_name, file_split_count, dependency_name):
        #split data into mRNA and rRNA so we can focus on the mRNA for the remainder of the analysis steps
        dep_loc = os.getcwd() + "/" + dependency_name + "/" + "data/final_results/"
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"

        orphan_split_folder = data_folder + "orphans/orphans_fastq/"
        pair_1_split_folder = data_folder + "pair_1/pair_1_fastq/"
        pair_2_split_folder = data_folder + "pair_2/pair_2_fastq/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(orphan_split_folder)
        self.make_folder(pair_1_split_folder)
        self.make_folder(pair_2_split_folder)

        file_splitter_orphans = self.tool_path_obj.Python + " " + self.tool_path_obj.File_splitter + " "
        file_splitter_orphans += dep_loc + "orphans_no_vectors.fastq "
        file_splitter_orphans += orphan_split_folder + "orphans "
        file_splitter_orphans += str(file_split_count)

        file_splitter_pair_1 = self.tool_path_obj.Python + " " + self.tool_path_obj.File_splitter + " "
        file_splitter_pair_1 += dep_loc + "pair_1_no_vectors.fastq "
        file_splitter_pair_1 += pair_1_split_folder + "pair_1 "
        file_splitter_pair_1 += str(file_split_count)

        file_splitter_pair_2 = self.tool_path_obj.Python + " " + self.tool_path_obj.File_splitter + " "
        file_splitter_pair_2 += dep_loc + "pair_2_no_vectors.fastq "
        file_splitter_pair_2 += pair_2_split_folder + "pair_2 "
        file_splitter_pair_2 += str(file_split_count)


        COMMANDS_rRNA_prep = [
            file_splitter_orphans,
            file_splitter_pair_1,
            file_splitter_pair_2
        ]

        return COMMANDS_rRNA_prep

    def create_rRNA_filter_command(self, stage_name, category, segment_root_name):
        #converts the fastq segments to fasta for infernal,
        #then takes the fasta segments, filters out the rRNA
        #then merges the split fastqs back together
        #called by each split file
        #category -> orphans, pair 1, pair 2
        #stage_name -> "rRNA_Filter"
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/" + category + "/"
        infernal_out_folder = data_folder + category + "_infernal/"
        mRNA_folder = data_folder + category + "_mRNA/"
        rRNA_folder = data_folder + category + "_rRNA/"
        fasta_folder = data_folder + category + "_fasta/"
        fastq_folder = data_folder + category + "_fastq/"

        infernal_out = infernal_out_folder + segment_root_name + ".infernal_out"
        fastq_in = fastq_folder + segment_root_name + ".fastq"
        fasta_io = fasta_folder + segment_root_name + ".fasta"

        self.make_folder(fasta_folder)
        self.make_folder(infernal_out_folder)
        self.make_folder(mRNA_folder)
        self.make_folder(rRNA_folder)

        convert_fastq_to_fasta = ">&2 echo converting " + category + " split file to fasta | "
        convert_fastq_to_fasta += self.tool_path_obj.vsearch
        convert_fastq_to_fasta += " --fastq_filter " + fastq_in
        convert_fastq_to_fasta += " --fastq_ascii " + self.Qual_str
        convert_fastq_to_fasta += " --fastaout " + fasta_io
        #print("converting", fastq_in )
        #print("placing in", fasta_io )

        infernal_command = ">&2 echo running infernal on " + category + " split file | "
        infernal_command += self.tool_path_obj.Infernal
        infernal_command += " -o /dev/null --tblout "
        infernal_command += infernal_out
        infernal_command += " --anytrunc --rfam -E 0.001 "
        infernal_command += self.tool_path_obj.Rfam + " "
        infernal_command += fasta_io

        rRNA_filtration = self.tool_path_obj.Python + " "
        rRNA_filtration += self.tool_path_obj.rRNA_filter + " "
        rRNA_filtration += infernal_out + " "
        rRNA_filtration += fastq_in + " "
        rRNA_filtration += mRNA_folder + " "
        rRNA_filtration += rRNA_folder

        COMMANDS_infernal = [
            convert_fastq_to_fasta,
            infernal_command,
            rRNA_filtration
        ]
        return COMMANDS_infernal

    def create_rRNA_filter_post_command(self, stage_name):
        #rRNA filtration orphaned some reads in the pairs.  We need to refilter the orphans.
        #Cat then refilter
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        pre_filter_folder = data_folder + "0_pre_orphans/"
        pre_filter_mRNA_folder = pre_filter_folder + "mRNA/"
        pre_filter_rRNA_folder = pre_filter_folder + "rRNA/"
        final_folder = data_folder + "final_results/"
        final_mRNA_folder = final_folder + "mRNA/"
        final_rRNA_folder = final_folder + "rRNA/"
        orphans_mRNA_folder = data_folder + "orphans/orphans_mRNA/"
        orphans_rRNA_folder = data_folder + "orphans/orphans_rRNA/"
        pair_1_mRNA_folder = data_folder + "pair_1/pair_1_mRNA/"
        pair_1_rRNA_folder = data_folder + "pair_1/pair_1_rRNA/"
        pair_2_mRNA_folder = data_folder + "pair_2/pair_2_mRNA/"
        pair_2_rRNA_folder = data_folder + "pair_2/pair_2_rRNA/"

        self.make_folder(pre_filter_folder)
        self.make_folder(pre_filter_mRNA_folder)
        self.make_folder(pre_filter_rRNA_folder)
        self.make_folder(final_folder)
        self.make_folder(final_mRNA_folder)
        self.make_folder(final_rRNA_folder)

        cat_orphans_mRNA = "cat " + orphans_mRNA_folder + "* 1>>" + pre_filter_mRNA_folder + "orphans.fastq"
        cat_orphans_rRNA = "cat " + orphans_rRNA_folder + "* 1>>" + pre_filter_rRNA_folder + "orphans.fastq"

        cat_pair_1_mRNA = "cat " + pair_1_mRNA_folder + "* 1>>" + pre_filter_mRNA_folder + "pair_1.fastq"
        cat_pair_1_rRNA = "cat " + pair_1_rRNA_folder + "* 1>>" + pre_filter_rRNA_folder + "pair_1.fastq"

        cat_pair_2_mRNA = "cat " + pair_2_mRNA_folder + "* 1>>" + pre_filter_mRNA_folder + "pair_2.fastq"
        cat_pair_2_rRNA = "cat " + pair_2_rRNA_folder + "* 1>>" + pre_filter_rRNA_folder + "pair_2.fastq"

        orphan_mRNA_filter = ">&2 echo filtering mRNA for orphans | "
        orphan_mRNA_filter += self.tool_path_obj.Python + " "
        orphan_mRNA_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_mRNA_filter += pre_filter_mRNA_folder + "pair_1.fastq "
        orphan_mRNA_filter += pre_filter_mRNA_folder + "pair_2.fastq "
        orphan_mRNA_filter += pre_filter_mRNA_folder + "orphans.fastq "
        orphan_mRNA_filter += final_mRNA_folder + "pair_1.fastq "
        orphan_mRNA_filter += final_mRNA_folder + "pair_2.fastq "
        orphan_mRNA_filter += final_mRNA_folder + "orphans.fastq"

        orphan_rRNA_filter = ">&2 echo filtering rRNA for orphans | "
        orphan_rRNA_filter += self.tool_path_obj.Python + " "
        orphan_rRNA_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_rRNA_filter += pre_filter_rRNA_folder + "pair_1.fastq "
        orphan_rRNA_filter += pre_filter_rRNA_folder + "pair_2.fastq "
        orphan_rRNA_filter += pre_filter_rRNA_folder + "orphans.fastq "
        orphan_rRNA_filter += final_rRNA_folder + "pair_1.fastq "
        orphan_rRNA_filter += final_rRNA_folder + "pair_2.fastq "
        orphan_rRNA_filter += final_rRNA_folder + "orphans.fastq"


        COMMANDS_rRNA_post = [
            cat_orphans_mRNA,
            cat_orphans_rRNA,
            cat_pair_1_mRNA,
            cat_pair_1_rRNA,
            cat_pair_2_mRNA,
            cat_pair_2_rRNA,
            orphan_mRNA_filter,
            orphan_rRNA_filter
        ]

        return COMMANDS_rRNA_post

    def create_repop_command(self, stage_name, preprocess_stage_name, dependency_stage_name):
        #This stage reintroduces the duplicate reads into the data.  We need it to count towards things.
        #Due to time, and hierarchical importance, we're leaving this stage alone.
        #Leaving it alone in a tangled state
        #But the issue is that by leaving it alone, we violate the design plan
        #The fix? We have to detect if preprocess has been run.  If so, pull the missing data there
        #if not,
        #What has to happen here:
        #-> detect if we've run the preprocess stage.
        #-> if it's run, grab data
        #-> if not, run our own custom preprocess up to what we need
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/" + "data/final_results/"
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        repop_folder = data_folder + "0_repop/"
        final_folder = data_folder + "final_results/"
        preprocess_subfolder = os.getcwd() + "/" + preprocess_stage_name + "/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(repop_folder)
        self.make_folder(final_folder)


        #we ran a previous preprocess.  grab files
        #need 3, 5(clstr only), and mRNA from the 2nd stage.
        #for the mRNA, we don't really care if it is.  This stage is just supposed to add in the missing duplicates from something that was stripped.

        hq_path = preprocess_subfolder + "data/4_orphan_read_filter/"
        cluster_path = preprocess_subfolder + "data/5_remove_duplicates/"

        repop_orphans = ">&2 echo Duplication repopulate Orphans | "
        repop_orphans += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_orphans += hq_path + "orphans.fastq" + " "   #in -> way back when things were quality-filtered.
                                                           #      step 2 in preprocess.  could still contain rRNA
        repop_orphans += dep_loc + "mRNA/orphans.fastq" + " "      #in -> rRNA filtration output
        repop_orphans += cluster_path + "orphans_unique.fastq.clstr" + " " #in -> duplicates filter output
        repop_orphans += repop_folder + "orphans.fastq"        #out


        repop_pair_1 = ">&2 echo Duplication repopulate pair 1 | "
        repop_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_1 += hq_path + "pair_1_match.fastq" + " "
        repop_pair_1 += dep_loc + "mRNA/pair_1.fastq" + " "
        repop_pair_1 += cluster_path + "pair_1_unique.fastq.clstr" + " "
        repop_pair_1 += repop_folder + "pair_1.fastq"

        repop_pair_2 = ">&2 echo Duplication repopulate pair 2 | "
        repop_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.duplicate_repopulate + " "
        repop_pair_2 += hq_path + "pair_2_match.fastq" + " "
        repop_pair_2 += dep_loc + "mRNA/pair_2.fastq" + " "
        repop_pair_2 += cluster_path + "pair_2_unique.fastq.clstr" + " "
        repop_pair_2 += repop_folder + "pair_2.fastq"

        orphan_repop_filter = ">&2 echo filtering mRNA for orphans | "
        orphan_repop_filter += self.tool_path_obj.Python + " "
        orphan_repop_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_repop_filter += repop_folder + "pair_1.fastq "
        orphan_repop_filter += repop_folder + "pair_2.fastq "
        orphan_repop_filter += repop_folder + "orphans.fastq "
        orphan_repop_filter += final_folder + "pair_1.fastq "
        orphan_repop_filter += final_folder + "pair_2.fastq "
        orphan_repop_filter += final_folder + "orphans.fastq"


        COMMANDS_Combine = [
        repop_orphans,
        repop_pair_1,
        repop_pair_2,
        orphan_repop_filter

        ]
        return COMMANDS_Combine

    def create_assemble_contigs_command(self, stage_name, dependency_stage_name):
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        spades_folder = data_folder + "0_spades/"
        bwa_folder = data_folder + "1_bwa_align/"
        sam_trimmer_folder = data_folder + "2_clean_sam/"
        mapped_reads_folder = data_folder + "3_mapped_reads/"
        final_folder = data_folder + "final_results/"
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(spades_folder)
        self.make_folder(bwa_folder)
        self.make_folder(sam_trimmer_folder)
        self.make_folder(mapped_reads_folder)
        self.make_folder(final_folder)

        #this assembles contigs
        spades = ">&2 echo Spades Contig assembly | "
        spades += self.tool_path_obj.Python + " "
        spades += self.tool_path_obj.Spades + " --rna"
        spades += " -1 " + dep_loc + "pair_1.fastq"  # in1 (pair 1)
        spades += " -2 " + dep_loc + "pair_2.fastq"  # in2 (pair 2)
        spades += " -s " + dep_loc + "orphans.fastq"  # in_single (orphans)
        spades += " -o " + spades_folder  # out

        spades_rename = "cp " + spades_folder + "transcripts.fasta" + " " + spades_folder + "contigs.fasta" # rename output

        bwa_index = self.tool_path_obj.BWA + " index -a bwtsw " + spades_folder + "contigs.fasta"

        #calls BWA, then uses SAMTools to get a report
        bwa_pair_1_contigs = ">&2 echo BWA pair contigs | "
        bwa_pair_1_contigs += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_pair_1_contigs += spades_folder + "contigs.fasta" + " "
        bwa_pair_1_contigs += dep_loc + "pair_1.fastq"
        bwa_pair_1_contigs += " > " + bwa_folder + "pair_1.sam"

        bwa_pair_2_contigs = ">&2 echo BWA pair contigs | "
        bwa_pair_2_contigs += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_pair_2_contigs += spades_folder + "contigs.fasta" + " "
        bwa_pair_2_contigs += dep_loc + "pair_2.fastq"
        bwa_pair_2_contigs += " > " + bwa_folder + "pair_2.sam"

        bwa_orphans_contigs = ">&2 echo BWA orphan contigs | "
        bwa_orphans_contigs += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 "
        bwa_orphans_contigs += spades_folder + "contigs.fasta" + " "
        bwa_orphans_contigs += dep_loc + "orphans.fastq"
        bwa_orphans_contigs += " > " + bwa_folder + "orphans.sam"

        sam_trimmer_orphans = ">&2 echo cleaning up orphans sam | "
        sam_trimmer_orphans += self.tool_path_obj.Python + " " + self.tool_path_obj.sam_trimmer + " "
        sam_trimmer_orphans += bwa_folder + "orphans.sam "
        sam_trimmer_orphans += sam_trimmer_folder + "orphans.sam"

        sam_trimmer_pair_1 = ">&2 echo cleaning up pair 1 sam | "
        sam_trimmer_pair_1 += self.tool_path_obj.Python + " " + self.tool_path_obj.sam_trimmer + " "
        sam_trimmer_pair_1 += bwa_folder + "pair_1.sam "
        sam_trimmer_pair_1 += sam_trimmer_folder + "pair_1.sam"

        sam_trimmer_pair_2 = ">&2 echo cleaning up pair 2 sam | "
        sam_trimmer_pair_2 += self.tool_path_obj.Python + " " + self.tool_path_obj.sam_trimmer + " "
        sam_trimmer_pair_2 += bwa_folder + "pair_2.sam "
        sam_trimmer_pair_2 += sam_trimmer_folder + "pair_2.sam"

        contig_duplicate_remover = ">&2 echo Removing consumed contigs from data | "
        contig_duplicate_remover += self.tool_path_obj.Python + " " + self.tool_path_obj.contig_duplicate_remover + " "
        contig_duplicate_remover += dep_loc + "pair_1.fastq "
        contig_duplicate_remover += dep_loc + "pair_2.fastq "
        contig_duplicate_remover += dep_loc + "orphans.fastq "
        contig_duplicate_remover += sam_trimmer_folder + "pair_1.sam "
        contig_duplicate_remover += sam_trimmer_folder + "pair_2.sam "
        contig_duplicate_remover += sam_trimmer_folder + "orphans.sam "
        contig_duplicate_remover += mapped_reads_folder

        map_read_contig_v2 = ">&2 echo map read contig v2 | "
        map_read_contig_v2 += self.tool_path_obj.Python + " " + self.tool_path_obj.map_read_contig_v2 + " "
        map_read_contig_v2 += dep_loc + "pair_1.fastq "
        map_read_contig_v2 += dep_loc + "pair_2.fastq "
        map_read_contig_v2 += dep_loc + "orphans.fastq "
        map_read_contig_v2 += sam_trimmer_folder + "pair_1.sam "
        map_read_contig_v2 += sam_trimmer_folder + "pair_2.sam "
        map_read_contig_v2 += sam_trimmer_folder + "orphans.sam "
        map_read_contig_v2 += final_folder + "contig_map.tsv"

        move_contigs = ">&2 echo Moving contigs to final folder | "
        move_contigs += "cp " + spades_folder + "contigs.fasta " + final_folder

        orphan_assembly_filter = ">&2 echo filtering paired reads for orphans | "
        orphan_assembly_filter += self.tool_path_obj.Python + " "
        orphan_assembly_filter += self.tool_path_obj.orphaned_read_filter + " "
        orphan_assembly_filter += mapped_reads_folder + "pair_1.fastq "
        orphan_assembly_filter += mapped_reads_folder + "pair_2.fastq "
        orphan_assembly_filter += mapped_reads_folder + "orphans.fastq "
        orphan_assembly_filter += final_folder + "pair_1.fastq "
        orphan_assembly_filter += final_folder + "pair_2.fastq "
        orphan_assembly_filter += final_folder + "orphans.fastq"

        sort_paired = ">&2 echo sorting paired reads | "
        sort_paired += self.tool_path_obj.Python + " " + self.tool_path_obj.sort_reads + " "
        sort_paired += final_folder + "pair_1.fastq" + " "
        sort_paired += final_folder + "pair_1_sorted.fastq" + " | "
        sort_paired += self.tool_path_obj.Python + " " + self.tool_path_obj.sort_reads + " "
        sort_paired += final_folder + "pair_2.fastq" + " "
        sort_paired += final_folder + "pair_2_sorted.fastq"

        """
        contig_merge = ">&2 echo Contig merge | "
        contig_merge += self.tool_path_obj.Python + " " + self.tool_path_obj.Map_reads_contigs + " " 
        contig_merge += self.Input_File1 + "_all_mRNA.fastq" + " " 
        contig_merge += self.Input_File2 + "_all_mRNA.fastq" + " " 
        contig_merge += self.Input_Filepath + "_all_mRNA_unpaired.fastq" + " " 
        contig_merge += self.Input_Filepath + "_contig_paired.sam" + " " 
        contig_merge += self.Input_Filepath + "_contig_unpaired.sam" + " " 
        contig_merge += self.Input_Filepath + "_contig_map.tsv"
        """
        COMMANDS_Assemble = [
                        #"mkdir -p " + os.path.join(self.Input_Path, os.path.splitext(self.Input_FName)[0]) + "_SpadesOut",
                        spades,
                        spades_rename,
                        bwa_index,
                        bwa_pair_1_contigs,
                        bwa_pair_2_contigs,
                        bwa_orphans_contigs,
                        sam_trimmer_orphans,
                        sam_trimmer_pair_1,
                        sam_trimmer_pair_2,
                        contig_duplicate_remover,
                        map_read_contig_v2,
                        move_contigs,
                        #bwa_unpaired_contigs#,
                        #contig_merge
                        orphan_assembly_filter,
                        sort_paired
                        ]
        return COMMANDS_Assemble
        
    def create_BWA_annotate_command(self, stage_name, dependency_stage_name, section):
        # meant to be called multiple times: section -> contigs, orphans, pair_1, pair_2
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        bwa_folder = data_folder + "0_bwa/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)

        if section == "contigs":
            section_file = section + ".fasta"
        else:
            section_file = section + ".fastq"
        bwa_job = ">&2 echo BWA on " + section + " | "
        bwa_job += self.tool_path_obj.BWA + " mem -t " + self.Threads_str + " " 
        bwa_job += self.tool_path_obj.DNA_DB + " " 
        bwa_job += dep_loc + section_file + " | " 
        bwa_job += self.tool_path_obj.SAMTOOLS + " view > " + bwa_folder + section + ".sam"
        
        
        
        return [bwa_job]
        
    def create_BWA_pp_command(self, stage_name, dependency_stage_name):    
        
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        bwa_folder = data_folder + "0_bwa/"
        final_folder = data_folder + "final_results/"
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(bwa_folder)
        self.make_folder(final_folder)
        
        map_read_bwa = ">&2 echo map read bwa v2 | "
        map_read_bwa += self.tool_path_obj.Python + " " + self.tool_path_obj.Map_reads_gene_BWA + " "
        map_read_bwa += self.tool_path_obj.DNA_DB + " "         #IN
        map_read_bwa += dep_loc + "contig_map.tsv "             #IN
        map_read_bwa += final_folder + "gene_map.tsv "          #OUT
        map_read_bwa += dep_loc + "contigs.fasta "              #IN
        map_read_bwa += bwa_folder + "contigs.sam "             #IN
        map_read_bwa += final_folder + "contigs.fasta" + " "    #OUT
        map_read_bwa += dep_loc + "orphans.fastq" + " "         #IN
        map_read_bwa += bwa_folder + "orphans.sam" + " "        #IN
        map_read_bwa += final_folder + "orphans.fasta" + " "    #OUT
        map_read_bwa += dep_loc + "pair_1.fastq" + " "          #IN
        map_read_bwa += bwa_folder + "pair_1.sam" + " "         #IN
        map_read_bwa += final_folder + "pair_1.fasta" + " "     #OUT
        map_read_bwa += dep_loc + "pair_2.fastq" + " "          #IN
        map_read_bwa += bwa_folder + "pair_2.sam" + " "         #IN
        map_read_bwa += final_folder + "pair_2.fasta"           #OUT
        
        move_contig_map = ">&2 echo copy contig map | "
        move_contig_map += "cp " + dep_loc + "contig_map.tsv " + final_folder + "contig_map.tsv"


        COMMANDS_Annotate_BWA = [
            map_read_bwa,
            move_contig_map
        ]
        return COMMANDS_Annotate_BWA

    def create_BLAT_annotate_command(self, stage_name, dependency_stage_name, section, fasta):
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        blat_folder = data_folder + "0_blat/"
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_folder)

        #blat_command = ">&2 echo BLAT annotation for " + section + " " + fasta +" | "
        blat_command = self.tool_path_obj.BLAT + " -noHead -minIdentity=90 -minScore=65 "
        blat_command += self.tool_path_obj.DNA_DB_Split + fasta + " "
        blat_command += dep_loc + section + ".fasta"
        blat_command += " -fine -q=rna -t=dna -out=blast8 -threads=2" + " "
        blat_command += blat_folder + section + "_" + fasta + ".blatout"

        return [blat_command]

    def create_BLAT_cat_command(self, stage_name, section):
        # this is meant to be called for each section: contigs, orphans, pair_1, pair_2
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        blat_folder = data_folder + "0_blat/"
        blat_merge_folder = data_folder + "1_blat_merge/"
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(blat_merge_folder)

        cat_command = "cat " + blat_folder + section + "*.blatout" + " > " + blat_merge_folder + section + ".blatout"
        return [cat_command]

    def create_BLAT_pp_command(self, stage_name, dependency_stage_name):
        # this call is meant to be run after the BLAT calls have been completed.
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/" # implied to be BWA
        blat_merge_folder = data_folder + "1_blat_merge/"
        
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)

        blat_pp = ">&2 echo BLAT post-processing | "
        blat_pp += self.tool_path_obj.Python + " " + self.tool_path_obj.Map_reads_gene_BLAT + " "
        blat_pp += self.tool_path_obj.DNA_DB + " "
        blat_pp += dep_loc + "contig_map.tsv" + " "
        blat_pp += dep_loc + "gene_map.tsv" + " "
        blat_pp += final_folder + "genes.fna" + " "
        blat_pp += final_folder + "gene_map.tsv "
        blat_pp += final_folder + "genes.fna "
        blat_pp += dep_loc + "contigs.fasta" + " "
        blat_pp += blat_merge_folder + "contigs.blatout" + " "
        blat_pp += final_folder + "contigs.fasta" + " "
        blat_pp += dep_loc + "orphans.fasta" + " "
        blat_pp += blat_merge_folder + "orphans.blatout" + " "
        blat_pp += final_folder + "orphans.fasta" + " "
        blat_pp += dep_loc + "pair_1.fasta" + " "
        blat_pp += blat_merge_folder + "pair_1.blatout" + " "
        blat_pp += final_folder + "pair_1.fasta" + " "
        blat_pp += dep_loc + "pair_2.fasta" + " "
        blat_pp += blat_merge_folder + "pair_2.blatout" + " "
        blat_pp += final_folder + "pair_2.fasta"
        
        move_contig_map = ">&2 echo copy contig map | "
        move_contig_map += "cp " + dep_loc + "contig_map.tsv " + final_folder + "contig_map.tsv"


        COMMANDS_Annotate_BLAT_Post = [
                        blat_pp,
                        move_contig_map
                        ]
        return COMMANDS_Annotate_BLAT_Post

    def create_DIAMOND_annotate_command(self, stage_name, dependency_stage_name, section):
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        diamond_folder = data_folder + "0_diamond/"
        section_folder = data_folder + section + "/"
        section_temp_folder = section_folder + "temp" + "/"
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(diamond_folder)
        self.make_folder(section_folder)
        self.make_folder(section_temp_folder)
        
        #COMMANDS_Annotate_Diamond = []
        #for j in range(0, 4):
        #    for i in range(1, count+1):
        #tag = "_dmnd_tmp_" + str(count)
        
        #Diamond_command_list = []
        #"mkdir -p " + os.path.splitext(self.Input_FName)[0] + tag,
        diamond_annotate = ">&2 echo gene annotate DIAMOND " + section +  " | "
        diamond_annotate += self.tool_path_obj.DIAMOND 
        diamond_annotate += " blastx -p " + self.Threads_str 
        diamond_annotate += " -d " + self.tool_path_obj.Prot_DB 
        diamond_annotate += " -q " + dep_loc + section + ".fasta" 
        diamond_annotate += " -o " +  diamond_folder + section + ".dmdout" 
        diamond_annotate += " -f 6 -t " + section_temp_folder
        diamond_annotate += " -k 10 --id 85 --query-cover 65 --min-score 60"            
                
         
        return [diamond_annotate]
        
    def create_DIAMOND_pp_command(self, stage_name, dependency_0_stage_name):    
        # the command just calls the merger program

        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc_0 = os.getcwd() + "/" + dependency_0_stage_name + "/data/final_results/" #implied to be blat pp
        diamond_folder = data_folder + "0_diamond/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)

        diamond_pp = ">&2 echo DIAMOND post process | "
        diamond_pp += self.tool_path_obj.Python         + " " + self.tool_path_obj.Map_reads_prot_DMND + " "
        diamond_pp += self.tool_path_obj.Prot_DB_plain  + " "                                                   #IN
        diamond_pp += dep_loc_0                         + "contig_map.tsv" + " "                                #IN
        diamond_pp += dep_loc_0                         + "gene_map.tsv" + " "                                  #IN
        diamond_pp += final_folder                      + "gene_map.tsv" + " "                                  #OUT
        diamond_pp += dep_loc_0                         + "genes.fna" + " "                                     #IN
        diamond_pp += final_folder                      + "proteins.faa" + " "                                  #OUT
        diamond_pp += dep_loc_0                         + "contigs.fasta" + " "                                 #IN
        diamond_pp += diamond_folder                    + "contigs.dmdout" + " "                                #IN
        diamond_pp += final_folder                      + "contigs.fasta" + " "                                 #OUT
        diamond_pp += dep_loc_0                         + "orphans.fasta" + " "                                 #IN
        diamond_pp += diamond_folder                    + "orphans.dmdout" + " "                                #IN
        diamond_pp += final_folder                      + "orphans.fasta" + " "                                 #OUT
        diamond_pp += dep_loc_0                         + "pair_1.fasta" + " "                                  #IN
        diamond_pp += diamond_folder                    + "pair_1.dmdout" + " "                                 #IN
        diamond_pp += final_folder                      + "pair_1.fasta" + " "                                  #OUT
        diamond_pp += dep_loc_0                         + "pair_2.fasta" + " "                                  #IN
        diamond_pp += diamond_folder                    + "pair_2.dmdout" + " "                                 #IN
        diamond_pp += final_folder                      + "pair_2.fasta"                                        #OUT


        COMMANDS_Annotate_Diamond_Post = [
        diamond_pp
        ]
        return COMMANDS_Annotate_Diamond_Post

    def create_taxonomic_annotation_command(self, current_stage_name, assemble_contigs_stage, diamond_stage):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        assemble_contigs_folder = os.getcwd() + "/" + assemble_contigs_stage + "/data/final_results/"
        diamond_folder = os.getcwd() + "/" + diamond_stage + "/data/final_results/"
        ga_taxa_folder = data_folder + "0_gene_taxa/"
        kaiju_folder = data_folder + "1_kaiju/"
        centrifuge_folder = data_folder + "2_centrifuge/"
        wevote_folder = data_folder + "3_wevote/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(ga_taxa_folder)
        self.make_folder(kaiju_folder)
        self.make_folder(centrifuge_folder)
        self.make_folder(wevote_folder)
        self.make_folder(final_folder)

        get_taxa_from_gene = ">&2 echo get taxa from gene | "
        get_taxa_from_gene += self.tool_path_obj.Python + " " + self.tool_path_obj.Annotated_taxid + " " # SLOW STEP
        get_taxa_from_gene += diamond_folder + "gene_map.tsv" + " "
        get_taxa_from_gene += self.tool_path_obj.accession2taxid + " "
        get_taxa_from_gene += ga_taxa_folder + "ga_taxon.tsv"

        kaiju_on_contigs = ">&2 echo kaiju on contigs | "
        kaiju_on_contigs += self.tool_path_obj.Kaiju
        kaiju_on_contigs += " -t " + self.tool_path_obj.nodes
        kaiju_on_contigs += " -f " + self.tool_path_obj.Kaiju_db
        kaiju_on_contigs += " -i " + assemble_contigs_folder + "contigs.fasta"
        kaiju_on_contigs += " -z " + self.Threads_str
        kaiju_on_contigs += " -o " + kaiju_folder + "contigs.tsv"

        kaiju_on_orphans = ">&2 echo kaiju on orphans | "
        kaiju_on_orphans += self.tool_path_obj.Kaiju
        kaiju_on_orphans += " -t " + self.tool_path_obj.nodes
        kaiju_on_orphans += " -f " + self.tool_path_obj.Kaiju_db
        kaiju_on_orphans += " -i " + assemble_contigs_folder + "orphans.fastq"
        kaiju_on_orphans += " -z " + self.Threads_str
        kaiju_on_orphans += " -o " + kaiju_folder + "orphans.tsv"

        kaiju_on_paired = ">&2 echo kaiju on pairs | "
        kaiju_on_paired += self.tool_path_obj.Kaiju
        kaiju_on_paired += " -t " + self.tool_path_obj.nodes
        kaiju_on_paired += " -f " + self.tool_path_obj.Kaiju_db
        kaiju_on_paired += " -i " + assemble_contigs_folder + "pair_1.fastq"
        kaiju_on_paired += " -j " + assemble_contigs_folder + "pair_2.fastq"
        kaiju_on_paired += " -z " + self.Threads_str
        kaiju_on_paired += " -o " + kaiju_folder + "pairs.tsv"

        cat_kaiju = ">&2 echo merging all kaiju results | "
        cat_kaiju += "cat "
        cat_kaiju += kaiju_folder + "contigs.tsv" + " "
        cat_kaiju += kaiju_folder + "orphans.tsv" + " "
        cat_kaiju += kaiju_folder + "pairs.tsv"
        cat_kaiju += " > " + kaiju_folder + "merged_kaiju.tsv"

        centrifuge_on_orphans = ">&2 echo centrifuge on orphans | "
        centrifuge_on_orphans += self.tool_path_obj.Centrifuge
        centrifuge_on_orphans += " -x " + self.tool_path_obj.Centrifuge_db
        centrifuge_on_orphans += " -1 " + assemble_contigs_folder + "pair_1.fastq"
        centrifuge_on_orphans += " -2 " + assemble_contigs_folder + "pair_2.fastq"
        centrifuge_on_orphans += " -U " + assemble_contigs_folder + "orphans.fastq"
        centrifuge_on_orphans += " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID"
        centrifuge_on_orphans += " --phred" + self.Qual_str
        centrifuge_on_orphans += " -p 6" #+ str(int(self.Threads_str)/2)
        centrifuge_on_orphans += " -S " + centrifuge_folder + "reads.tsv"
        centrifuge_on_orphans += " --report-file " + centrifuge_folder + "reads.txt"

        centrifuge_on_contigs = ">&2 echo centrifuge on contigs | "
        centrifuge_on_contigs += self.tool_path_obj.Centrifuge
        centrifuge_on_contigs += " -f -x " + self.tool_path_obj.Centrifuge_db
        centrifuge_on_contigs += " -U " + assemble_contigs_folder + "contigs.fasta"
        centrifuge_on_contigs += " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID"
        centrifuge_on_contigs += " --phred" + self.Qual_str
        centrifuge_on_contigs += " -p 6" #+ str(int(self.Threads_str)/2)
        centrifuge_on_contigs += " -S " + centrifuge_folder + "contigs.tsv"
        centrifuge_on_contigs += " --report-file " + centrifuge_folder + "contigs.txt"
        
        cat_centrifuge = ">&2 echo combining all centrifuge results | "
        cat_centrifuge += "cat "
        cat_centrifuge += centrifuge_folder + "reads.tsv" + " "
        cat_centrifuge += centrifuge_folder + "contigs.tsv"
        cat_centrifuge += " > " + centrifuge_folder + "merged_centrifuge.tsv"

        wevote_combine = ">&2 echo combining classification outputs for wevote | "
        wevote_combine += self.tool_path_obj.Python + " " + self.tool_path_obj.Classification_combine + " "
        wevote_combine += assemble_contigs_folder + "contig_map.tsv"
        wevote_combine += " " + wevote_folder + "wevote_ensemble.csv" + " "
        wevote_combine += ga_taxa_folder + "ga_taxon.tsv" + " "
        wevote_combine += kaiju_folder + "merged_kaiju.tsv" + " "
        wevote_combine += centrifuge_folder + "merged_centrifuge.tsv"

        wevote_call = ">&2 echo Running WEVOTE | "
        wevote_call += self.tool_path_obj.WEVOTE
        wevote_call += " -i " + wevote_folder + "wevote_ensemble.csv"
        wevote_call += " -d " + self.tool_path_obj.WEVOTEDB
        wevote_call += " -p " + wevote_folder + "wevote"
        wevote_call += " -n " + self.Threads_str
        wevote_call += " -k " + "2"
        wevote_call += " -a " + "0"
        wevote_call += " -s " + "0"

        awk_cleanup = ">&2 echo AWK cleanup of WEVOTE results | "
        awk_cleanup += "awk -F \'\\t\' \'{print \"C\\t\"$1\"\\t\"$9}\' "
        awk_cleanup += wevote_folder + "wevote_WEVOTE_Details.txt"
        awk_cleanup += " > " + final_folder + "taxonomic_classifications.tsv"

        #taxid_to_english = ">&2 echo Constrain classification | "
        #taxid_to_english += self.tool_path_obj.Python + " " + self.tool_path_obj.Constrain_classification + " "
        #taxid_to_english += "family" + " "
        #taxid_to_english += self.Input_Filepath + "_WEVOTEOut.tsv" + " "
        #taxid_to_english += self.tool_path_obj.Nodes + " "
        #taxid_to_english += self.tool_path_obj.Names + " "
        #taxid_to_english += self.Input_Filepath + "_WEVOTEOut_family.tsv"

        #kaiju_to_krona = ">&2 echo kaiju to krona | "
        #kaiju_to_krona += self.tool_path_obj.Kaiju2krona
        #kaiju_to_krona += " -t " + self.tool_path_obj.nodes
        #kaiju_to_krona += " -n " + self.tool_path_obj.names
        #kaiju_to_krona += " -i " + self.Input_Filepath + "_WEVOTEOut_family.tsv"
        #kaiju_to_krona += " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt"
        
        #awk_cleanup_krona = "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        
        #kt_import_text_cleanup = self.tool_path_obj.ktImportText + " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"

        COMMANDS_Classify = [
            get_taxa_from_gene,
            kaiju_on_contigs,
            kaiju_on_orphans,
            kaiju_on_paired,
            cat_kaiju,
            centrifuge_on_orphans,
            centrifuge_on_contigs,
            cat_centrifuge,
            wevote_combine,
            wevote_call,
            awk_cleanup
            #taxid_to_english,
            #kaiju_to_krona,
            #awk_cleanup_krona,
            #kt_import_text_cleanup
            ]
        return COMMANDS_Classify 

    def create_EC_DETECT_prep(self, current_stage_name, diamond_stage, file_split_count):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        diamond_folder = os.getcwd() + "/" + diamond_stage + "/data/final_results/"
        proteins_folder = data_folder + "0_proteins/"
        detect_folder = data_folder + "1_detect/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(proteins_folder)
        self.make_folder(detect_folder)
        self.make_folder(final_folder)

        file_splitter = ">&2 echo splitting protein files for Detect | "
        file_splitter += self.tool_path_obj.Python + " " + self.tool_path_obj.File_splitter + " "
        file_splitter += diamond_folder + "proteins.faa" + " "
        file_splitter += proteins_folder + "protein" + " "
        file_splitter += str(file_split_count)

        COMMANDS_DETECT_prep = [
            file_splitter
            ]
        return COMMANDS_DETECT_prep

    def create_EC_DETECT_command(self, current_stage_name, prot_name):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/" + "/"
        proteins_folder = data_folder + "0_proteins/"
        detect_folder = data_folder + "1_detect/"
        prot_folder = detect_folder + prot_name + "/"

        self.make_folder(prot_folder)

        detect_protein = "(cd " + prot_folder + " && "
        detect_protein += ">&2 echo running detect on split file " + prot_name + " | "
        detect_protein += self.tool_path_obj.Python + " "
        detect_protein += self.tool_path_obj.Detect + " "
        detect_protein += proteins_folder + prot_name + ".fasta"
        detect_protein += " --output_file " + detect_folder + prot_name + ".detect"
        detect_protein += " --top_predictions_file " + detect_folder + prot_name + ".toppred"
        detect_protein += " --db " + self.tool_path_obj.DetectDB
        detect_protein += " --blastp " + self.tool_path_obj.Blastp
        detect_protein += " --needle " + self.tool_path_obj.Needle + ")"

        COMMANDS_DETECT = [
            detect_protein
            ]

        return COMMANDS_DETECT
        
    def create_EC_PRIAM_DIAMOND_command(self, current_stage_name, diamond_stage):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        diamond_folder = os.getcwd() + "/" + diamond_stage + "/data/final_results/"
        PRIAM_folder = data_folder + "2_priam/"
        diamond_ea_folder = data_folder + "3_diamond/"

        self.make_folder(PRIAM_folder)
        self.make_folder(diamond_ea_folder)

        PRIAM_command = ">&2 echo running PRIAM | "
        PRIAM_command += self.tool_path_obj.Java + " "
        PRIAM_command += self.tool_path_obj.Priam
        PRIAM_command += " -n " + "proteins_priam" + " "
        PRIAM_command += " -i " + diamond_folder + "proteins.faa"
        PRIAM_command += " -p " + self.tool_path_obj.PriamDB
        PRIAM_command += " -o " + PRIAM_folder
        PRIAM_command += " --np " + self.Threads_str
        PRIAM_command += " --bh --cc --cg --bp --bd "
        PRIAM_command += self.tool_path_obj.BLAST_dir

        diamond_ea_command = ">&2 echo running Diamond enzyme annotation | "
        diamond_ea_command += self.tool_path_obj.DIAMOND + " blastp"
        diamond_ea_command += " -p " + self.Threads_str
        diamond_ea_command += " --query " + diamond_folder + "proteins.faa"
        diamond_ea_command += " --db " + self.tool_path_obj.SWISS_PROT
        diamond_ea_command += " --outfmt " + "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident"
        diamond_ea_command += " --out " + diamond_ea_folder + "proteins.blastout"
        diamond_ea_command += " --evalue 0.0000000001 --max-target-seqs 1"

        COMMANDS_PRIAM_DIAMOND = [
            PRIAM_command,
            diamond_ea_command
            ]

        return COMMANDS_PRIAM_DIAMOND

    def create_EC_postprocess_command(self, current_stage_name, diamond_stage):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        diamond_folder = os.getcwd() + "/" + diamond_stage + "/data/final_results/"
        detect_folder = data_folder + "1_detect/"
        PRIAM_folder = data_folder + "2_priam/"
        diamond_ea_folder = data_folder + "3_diamond/"
        final_folder = data_folder + "final_results/"

        
        
        combine_detect = "cat " + detect_folder + "protein_*.toppred"
        combine_detect += " > " + detect_folder + "proteins.toppred"

        postprocess_command = ">&2 echo combining enzyme annotation output | "
        postprocess_command += self.tool_path_obj.Python + " "
        postprocess_command += self.tool_path_obj.EC_Annotation_Post + " "
        postprocess_command += diamond_folder + "proteins.faa" + " "
        postprocess_command += detect_folder + "proteins.toppred" + " "
        postprocess_command += os.path.join(PRIAM_folder, "PRIAM_proteins_priam", "ANNOTATION", "sequenceECs.txt") + " "
        postprocess_command += diamond_ea_folder + "proteins.blastout" + " "
        postprocess_command += self.tool_path_obj.SWISS_PROT + " "
        postprocess_command += self.tool_path_obj.SWISS_PROT_map + " "
        postprocess_command += final_folder


        COMMANDS_EC_Postprocess = [
            combine_detect,
            postprocess_command
            ]
                
        return COMMANDS_EC_Postprocess

    def create_Network_generation_command(self, current_stage_name, diamond_stage, taxonomic_annotation_stage, enzyme_annotation_stage):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        diamond_folder = os.getcwd() + "/" + diamond_stage + "/data/final_results/"
        ta_folder = os.getcwd() + "/" + taxonomic_annotation_stage + "/data/final_results/"
        ea_folder = os.getcwd() + "/" + enzyme_annotation_stage + "/data/final_results/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)

        network_generation = ">&2 echo generating RPKM table and Cytoscape network | "
        network_generation += self.tool_path_obj.Python + " "
        network_generation += self.tool_path_obj.RPKM + " "
        network_generation += self.tool_path_obj.nodes + " "
        network_generation += diamond_folder + "gene_map.tsv" + " "
        network_generation += ta_folder + "taxonomic_classifications.tsv" + " "
        network_generation += ea_folder + "proteins.ECs_All" + " "
        network_generation += final_folder + "RPKM_table.tsv" + " "
        network_generation += final_folder + "Cytoscape_network.tsv" + " "

        COMMANDS_Network = [
            network_generation
            ]
        return COMMANDS_Network

    def create_visualization_command(self, current_stage_name, network_stage):
        subfolder = os.getcwd() + "/" + current_stage_name + "/"
        data_folder = subfolder + "data/"
        mpl_folder = data_folder + "0_MPL/"
        network_folder = os.getcwd() + "/" + network_stage + "/data/final_results/"
        final_folder = data_folder + "final_results/"

        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(mpl_folder)
        self.make_folder(final_folder)

        chart_generation = ">&2 echo generating visualizations for ECs | "
        chart_generation += self.tool_path_obj.Python + " "
        chart_generation += self.tool_path_obj.chart + " "
        chart_generation += network_folder + "Cytoscape_network.tsv" + " "
        chart_generation += mpl_folder

        final_chart = "cp " + mpl_folder + "All_EC.png"  + " " + final_folder + "All_EC.png"

        COMMANDS_visualization = [
            chart_generation,
            final_chart
            ]
        return COMMANDS_visualization