#The functions here generate the pipeline commands.
#Each command module is made up of (often) many sub stages that are used to get the final result.
#If you want to move around the ordering, you'd do that here.

import os
import sys
import mt_pipe_paths as mpp
import subprocess as sp
#------------------------------------------------------
# pragmas needed for command construction

PBS_Submit_LowMem = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=00:15:00
#PBS -N NAME
#PBS -e ERROR
#PBS -o OUTPUT

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras anaconda3/4.0.0
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

PBS_Submit_HighMem = """#!/bin/bash
#PBS -l nodes=1:m64g:ppn=16,walltime=12:00:00 -q sandy
#PBS -N NAME
#PBS -e ERROR
#PBS -o OUTPUT

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras anaconda3/4.0.0
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=16
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

PBS_Submit_vHighMem = """#!/bin/bash
#PBS -l nodes=1:m128g:ppn=16,walltime=1:00:00 -q sandy
#PBS -N NAME
#PBS -e ERROR
#PBS -o OUTPUT

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras anaconda3/4.0.0
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=20
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

#--------------------------------------------------------------------------
# class definition

class mt_pipe_commands:
    #--------------------------------------------------------------------
    # constructor:
    # there should only be one of these objects used for an entire pipeline.
    def __init__(self, Quality_score = 33, Thread_count = 16, raw_sequence_path_0 = None, raw_sequence_path_1 = None):
        # path to the raw genome sequence file
        Input_File = os.getcwd()
        if not(raw_sequence_path_0 is None):
            self.raw_sequence_path_0 = raw_sequence_path_0
            print("raw sequence 0:", self.raw_sequence_path_0)

        if not(raw_sequence_path_1 is None):
            self.raw_sequence_path_1 = raw_sequence_path_1
            print("raw seqeunce 1:", self.raw_sequence_path_1)
        
            
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.Input_File1 = self.Input_Filepath + "1"
        self.Input_File2 = self.Input_Filepath + "2"
        self.Input_FName = os.path.basename(Input_File)
        
        self.Qual_str = str(Quality_score)
        self.Input_Path = os.path.dirname(Input_File)
        self.Contigs = os.path.join(self.Input_Path, os.path.splitext(self.Input_FName)[0] + "_SpadesOut", "contigs.fasta")
        self.Threads_str= str(Thread_count)
        self.EC_Split = os.path.join(self.Input_Filepath + "_EC_Annotation", "Split")
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        self.Host_Contaminants = self.Input_Filepath + "_host_contaminents_seq.fasta"
        self.Vector_Contaminants = self.Input_Filepath + "_vector_contaminants_seq.fasta"
        
        print("input filepath:", self.Input_Filepath)
        print("input file 1:", self.Input_File1)
        print("input file 2:", self.Input_File2)
        print("input FName:", self.Input_FName)
        print("input path:", self.Input_Path)
        print("where are we now:", os.getcwd())
    
    #-----------------------------------------------------------
    # support functions
    #def make_stage_dir(self, stage_name):
    #    if not(os.path.exists(self.)):
    def make_folder(self, folder_path):
        if not(os.path.exists(folder_path)):
            os.makedirs(folder_path)
    
    def create_pbs_and_launch(self, job_name, command_list, mode = "low", dependency_list = None, run_job = False,  inner_name = None):
        #create the pbs job, and launch items
        #job name: string tag for export file name
        #command list:  list of command statements for writing
        #mode: selection of which pbs template to use: default -> low memory
        #dependency_list: if not empty, will append wait args to qsub subprocess call. it's polymorphic
        #returns back the job ID given from qsub
        
        pbs_template = ""
        if(mode == "med"):
            pbs_template = PBS_Submit_HighMem
        elif(mode == "high"):
            pbs_template = PBS_Submit_vHighMem
        else:
            pbs_template = PBS_Submit_LowMem
        
        
        pbs_script_full_path = os.getcwd() + "/" + job_name +"/" + job_name
        if(not inner_name is None):
            pbs_script_full_path = os.getcwd() + "/" + job_name + "/" + inner_name
            
        try:
            with open(pbs_script_full_path + ".pbs", "w+") as PBS_script_out:
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
                if (isinstance(dependency_list, int)):
                    #single dep
                    dep_str = "-W depend=afterok:" + str(dependency_list)
                    print(job_name, "running with single dependency")
                    
                elif(isinstance(dependency_list, list)):
                    # multiple deps
                    dep_str = "-W depend=afterok:"
                    for item in dependency_list:
                        dep_str += ":" + str(item)
                    print(job_name, "running with multiple dependencies")    
                        
                elif(dependency_list is None):
                    print(job_name, "running without dependency")
                else:
                    print("This isn't supposed to happen")
                    
                if(run_job): # a lock built for testing syntax, but not run
                    try:
                        if not dep_str == "":
                            print("dep string not empty")
                            job_id = sp.check_output(["qsub", pbs_script_full_path+".pbs", dep_str])
                        else:
                            job_id = sp.check_output(["qsub", pbs_script_full_path+".pbs"])
                        #return val is a binary string, so we convert, and extract only the numeric part
                        job_id = job_id.decode('ascii')
                        job_id = int(job_id.split('.')[0])
                    
                        return job_id
                    except Exception as e:
                        print("subprocess call error:", e)
                else:
                    return 0
                
                
        except Exception as e:
            # error catchall 
            print("Failure at pbs creation:", e)
            sys.exit()
    
    def create_pre_single_command(self, stage_name):
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        
        print("not ready")
        adapter_removal = mpp.AdapterRemoval + "--file1 " + ""
        COMMANDS_PRE = []
    
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
        sort_pair_1 += mpp.Python + " " + mpp.Sort_Reads + " "
        sort_pair_1 += self.raw_sequence_path_0
        sort_pair_1 += " " + sorted_raw_folder + "pair_1_sorted.fastq"
        
        sort_pair_2 = ">&2 echo Sorting pair 2 | "
        sort_pair_2 += mpp.Python + " " + mpp.Sort_Reads + " "
        sort_pair_2 += self.raw_sequence_path_1
        sort_pair_2 += " " + sorted_raw_folder + "pair_2_sorted.fastq"
        
        adapter_folder = data_folder + "1_adapter_removal/"
        self.make_folder(adapter_folder)
        adapter_removal_line = ">&2 echo Removing adapters | "
        adapter_removal_line += mpp.AdapterRemoval 
        #adapter_removal_line += " --file1 " + self.raw_sequence_path_0
        #adapter_removal_line += " --file2 " + self.raw_sequence_path_1
        adapter_removal_line += " --file1 " + sorted_raw_folder + "pair_1_sorted.fastq"
        adapter_removal_line += " --file2 " + sorted_raw_folder + "pair_2_sorted.fastq"
        adapter_removal_line += " --qualitybase " + str(self.Qual_str) #must be either 33 or 64
        adapter_removal_line += " --threads " + self.Threads_str  
        adapter_removal_line += " --minlength " + "30" 
        adapter_removal_line += " --basename " + adapter_folder + os.path.splitext(self.Input_FName)[0]  
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
        vsearch_merge += mpp.vsearch 
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
        vsearch_filter_0 += mpp.vsearch 
        vsearch_filter_0 += " --fastq_filter " + vsearch_merge_folder + "orphans.fastq" 
        vsearch_filter_0 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_0 += " --fastq_maxee " + "2.0" 
        vsearch_filter_0 += " --fastqout " + vsearch_filter_folder + "orphans_hq.fastq"
        
        #then move onto the standalones in pair 1
        vsearch_filter_1 = ">&2 echo low-quality filter on pair 1 | "
        vsearch_filter_1 += mpp.vsearch  
        vsearch_filter_1 += " --fastq_filter " + vsearch_merge_folder + "pair_1_merge_reject.fastq" 
        vsearch_filter_1 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_1 += " --fastq_maxee " + "2.0" 
        vsearch_filter_1 += " --fastqout " + vsearch_filter_folder + "pair_1_hq.fastq"
        
        vsearch_filter_2 = ">&2 echo low-quality filter on pair 2 | "
        vsearch_filter_2 += mpp.vsearch 
        vsearch_filter_2 += " --fastq_filter " + vsearch_merge_folder + "pair_2_merge_reject.fastq" 
        vsearch_filter_2 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_2 += " --fastq_maxee " + "2.0" 
        vsearch_filter_2 += " --fastqout " + vsearch_filter_folder + "pair_2_hq.fastq"
        
        #redistribute data into orphans, or paired-reads
        orphan_read_filter_folder = data_folder + "4_orphan_read_filter/"
        self.make_folder(orphan_read_filter_folder)
        orphan_read_filter = ">&2 echo moving newly orphaned hq reads | "
        orphan_read_filter += mpp.Python + " " 
        orphan_read_filter += mpp.orphaned_read_filter + " " 
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
        cdhit_orphans += mpp.cdhit_dup + " -i " 
        cdhit_orphans += orphan_read_filter_folder + "orphans.fastq" 
        cdhit_orphans += " -o " + cdhit_folder + "orphans_unique.fastq"
        
        
        #move_unpaired_cluster = "mv " 
        #move_unpaired_cluster += self.Input_Filepath + "_unpaired_unique.fastq.clstr " 
        #move_unpaired_cluster += self.Input_Filepath + "_unpaired.clstr"
        
        #remove duplicates in the pairs
        cdhit_pair_1 = ">&2 echo remove duplicates from pair 1 | "
        cdhit_pair_1 += mpp.cdhit_dup
        cdhit_pair_1 += " -i " + orphan_read_filter_folder + "pair_1_match.fastq"
        cdhit_pair_1 += " -o " + cdhit_folder + "pair_1_unique.fastq"
        
        cdhit_pair_2 = ">&2 echo remove duplicates from pair 2 | " 
        cdhit_pair_2 += mpp.cdhit_dup
        cdhit_pair_2 += " -i " + orphan_read_filter_folder + "pair_2_match.fastq"
        cdhit_pair_2 += " -o " + cdhit_folder + "pair_2_unique.fastq"
        
        #move_paired_cluster = "mv " + self.Input_File1 + "_paired_unique.fastq.clstr" + " " + self.Input_Filepath + "_paired.clstr"
        #NOTE: we overwrite host contaminants here. why? because the preprocess folder is created
        #in this function only
        host_removal_folder = data_folder + "6_remove_host/"
        self.make_folder(host_removal_folder)
        
        self.Host_Contaminants = host_removal_folder + "host_contaminents_seq.fasta"
        copy_host = ">&2 echo Copy the host file over | "
        copy_host += "cp " + mpp.Host + " " + self.Host_Contaminants
        
        #craft a BWA index for the host sequences
        bwa_hr_prep = ">&2 echo make host contaminants index for BWA | "
        bwa_hr_prep += mpp.BWA + " index -a bwtsw " + self.Host_Contaminants
        
        samtools_hr_prep = ">&2 echo SAMTOOLS host contaminant prep | "
        samtools_hr_prep += mpp.SAMTOOLS + " faidx " + self.Host_Contaminants
        
        #host removal on unique orphans
        bwa_hr_orphans = ">&2 echo BWA host remove on orphans | "
        bwa_hr_orphans += mpp.BWA + " mem -t " 
        bwa_hr_orphans += self.Threads_str + " " 
        bwa_hr_orphans += self.Host_Contaminants + " " 
        bwa_hr_orphans += cdhit_folder + "orphans_unique.fastq" 
        bwa_hr_orphans += " > " + host_removal_folder + "orphans_no_host.sam"
        
        #annoying type conversion pt 1
        samtools_hr_orphans_sam_to_bam = ">&2 echo convert orphans hr files pt1 | "
        samtools_hr_orphans_sam_to_bam += mpp.SAMTOOLS 
        samtools_hr_orphans_sam_to_bam += " view -bS " + host_removal_folder + "orphans_no_host.sam" 
        samtools_hr_orphans_sam_to_bam += " > " + host_removal_folder + "orphans_no_host.bam"
        #annoying type conversion pt 2
        samtools_no_host_orphans_bam_to_fastq = ">&2 echo convert orphans hr files pt2 | "
        samtools_no_host_orphans_bam_to_fastq += mpp.SAMTOOLS 
        samtools_no_host_orphans_bam_to_fastq += " fastq -n -f 4" + " -0 " + host_removal_folder + "orphans_no_host.fastq" + " "
        samtools_no_host_orphans_bam_to_fastq += host_removal_folder + "orphans_no_host.bam"
        
        #apparently, we're to keep the host separation
        samtools_host_orphans_bam_to_fastq = ">&2 echo convert orphans hr files pt3 | "
        samtools_host_orphans_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -F 4" 
        samtools_host_orphans_bam_to_fastq += " -0 " + host_removal_folder + "orphans_host_only.fastq" + " " 
        samtools_host_orphans_bam_to_fastq += host_removal_folder + "orphans_no_host.bam"
        
        #host-remove the rest
        bwa_hr_pair = ">&2 echo bwa pair host remove | "
        bwa_hr_pair += mpp.BWA + " mem -t " + self.Threads_str + " " + self.Host_Contaminants + " " 
        bwa_hr_pair += cdhit_folder + "pair_1_unique.fastq" + " " 
        bwa_hr_pair += cdhit_folder + "pair_2_unique.fastq" 
        bwa_hr_pair += " > " + host_removal_folder + "pair_no_host.sam"
        
        
        #separating bwa results back into paired reads
        samtools_host_pair_sam_to_bam = ">&2 echo convert pair hr files pt1 | "
        samtools_host_pair_sam_to_bam += mpp.SAMTOOLS + " view -bS " + host_removal_folder + "pair_no_host.sam"
        samtools_host_pair_sam_to_bam += " > " + host_removal_folder + "pair_no_host.bam"
        
        #stuff that doesn't match with the host
        samtools_no_host_pair_bam_to_fastq = ">&2 echo convert pair hr files pt2 | "
        samtools_no_host_pair_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -f 13" 
        samtools_no_host_pair_bam_to_fastq += " -1 " + host_removal_folder + "pair_1_no_host.fastq" # out
        samtools_no_host_pair_bam_to_fastq += " -2 " + host_removal_folder + "pair_2_no_host.fastq" # out
        samtools_no_host_pair_bam_to_fastq += " " + host_removal_folder + "pair_no_host.bam" #in
        
        #stuff that matches with the host (why keep it?  request from john)
        samtools_host_pair_bam_to_fastq = ">&2 echo convert pair hr files pt3 | "
        samtools_host_pair_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -F 4" 
        samtools_host_pair_bam_to_fastq += " -1 " + host_removal_folder + "pair_1_host_only.fastq" 
        samtools_host_pair_bam_to_fastq += " -2 " + host_removal_folder + "pair_2_host_only.fastq" 
        samtools_host_pair_bam_to_fastq += " " + host_removal_folder + "pair_no_host.bam"
        
        #blast prep
        make_blast_db_host = ">&2 echo Make BLAST db for host contaminants | "
        make_blast_db_host += mpp.Makeblastdb + " -in " + self.Host_Contaminants + " -dbtype nucl"

        vsearch_filter_3 = ">&2 echo Convert orphans for BLAT | "     
        vsearch_filter_3 += mpp.vsearch 
        vsearch_filter_3 += " --fastq_filter " + host_removal_folder + "orphans_no_host.fastq" 
        vsearch_filter_3 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_3 += " --fastaout " + host_removal_folder + "orphans_no_host.fasta"
        
        vsearch_filter_4 = ">&2 echo Convert pair 1 for BLAT | "
        vsearch_filter_4 += mpp.vsearch 
        vsearch_filter_4 += " --fastq_filter " + host_removal_folder + "pair_1_no_host.fastq" 
        vsearch_filter_4 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_4 += " --fastaout " + host_removal_folder + "pair_1_no_host.fasta"
        
        vsearch_filter_5 = ">&2 echo Convert pair 2 for BLAT | "
        vsearch_filter_5 += mpp.vsearch
        vsearch_filter_5 += " --fastq_filter " + host_removal_folder + "pair_2_no_host.fastq" 
        vsearch_filter_5 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_5 += " --fastaout " + host_removal_folder + "pair_2_no_host.fasta"
        
        blat_hr_orphans = ">&2 echo BLAT hr orphans | "
        blat_hr_orphans += mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " 
        blat_hr_orphans += self.Host_Contaminants + " " 
        blat_hr_orphans += host_removal_folder + "orphans_no_host.fasta" 
        blat_hr_orphans += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str 
        blat_hr_orphans += " " + host_removal_folder + "orphans_no_host.blatout"

        blat_hr_pair_1 = ">&2 echo BLAT hr pair 1 | "    
        blat_hr_pair_1 += mpp.BLAT 
        blat_hr_pair_1 += " -noHead -minIdentity=90 -minScore=65 " + self.Host_Contaminants + " " 
        blat_hr_pair_1 += host_removal_folder + "pair_1_no_host.fasta" 
        blat_hr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_1 += " " + host_removal_folder + "pair_1_no_host.blatout"
        
        blat_hr_pair_2 = ">&2 echo BLAT hr pair 2 | "
        blat_hr_pair_2 += mpp.BLAT 
        blat_hr_pair_2 += " -noHead -minIdentity=90 -minScore=65 " + self.Host_Contaminants + " " 
        blat_hr_pair_2 += host_removal_folder + "pair_2_no_host.fasta" 
        blat_hr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str
        blat_hr_pair_2 += " " + host_removal_folder + "pair_2_no_host.blatout"
        
        # HR BLAT
        blat_hr_folder = data_folder + "7_blat_hr/"
        self.make_folder(blat_hr_folder)
        hr_orphans = ">&2 echo BLAT contaminant orphans | "
        hr_orphans += mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " 
        hr_orphans += host_removal_folder + "orphans_no_host.fastq" + " " # in
        hr_orphans += host_removal_folder + "orphans_no_host.blatout" + " " #in 
        hr_orphans += blat_hr_folder + "orphans_no_host.fastq" + " " #out
        hr_orphans += blat_hr_folder + "orphans_host_only.fastq" #out
        
        hr_pair_1 = ">&2 echo BLAT contaminant pair 1 | "
        hr_pair_1 += mpp.Python + " " 
        hr_pair_1 += mpp.BLAT_Contaminant_Filter + " " 
        hr_pair_1 += host_removal_folder + "pair_1_no_host.fastq" + " " 
        hr_pair_1 += host_removal_folder + "pair_1_no_host.blatout" + " " 
        hr_pair_1 += blat_hr_folder + "pair_1_no_host.fastq" + " " 
        hr_pair_1 += blat_hr_folder + "pair_1_host_only.fastq"
        
        hr_pair_2 = ">&2 echo BLAT contaminant pair 2 | "
        hr_pair_2 += mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " 
        hr_pair_2 += host_removal_folder + "pair_2_no_host.fastq" + " " 
        hr_pair_2 += host_removal_folder + "pair_2_no_host.blatout" + " " 
        hr_pair_2 += blat_hr_folder + "pair_2_no_host.fastq" + " " 
        hr_pair_2 += blat_hr_folder + "pair_2_host_only.fastq"
        
        #-----------------------------
        # Vector removal
        # vectors are artifacts that came from the read-gathering
        vector_removal_folder = data_folder + "8_vector_removal/"
        self.make_folder(vector_removal_folder)
        self.Vector_Contaminants = vector_removal_folder + "vector_contaminants_seq.fasta"
        
        copy_vector = ">&2 echo copy vector prep | "
        copy_vector += "cp " + mpp.UniVec_Core + " " + self.Vector_Contaminants
        
        bwa_vr_prep = ">&2 echo BWA vector prep | "
        bwa_vr_prep += mpp.BWA + " index -a bwtsw " + self.Vector_Contaminants
        
        samtools_vr_prep = ">&2 echo samtools vector prep | "
        samtools_vr_prep += mpp.SAMTOOLS + " faidx " + self.Vector_Contaminants
        
        bwa_vr_orphans = ">&2 echo BWA vr oprhans | "
        bwa_vr_orphans += mpp.BWA + " mem -t " + self.Threads_str + " " 
        bwa_vr_orphans += self.Vector_Contaminants + " " 
        bwa_vr_orphans += blat_hr_folder + "orphans_no_host.fastq"  
        bwa_vr_orphans += " > " + vector_removal_folder + "orphans_no_vectors.sam"
        
        samtools_no_vector_orphans_sam_to_bam = ">&2 echo samtools vr oprhans pt 1 | "
        samtools_no_vector_orphans_sam_to_bam += mpp.SAMTOOLS + " view -bS " 
        samtools_no_vector_orphans_sam_to_bam += vector_removal_folder + "orphans_no_vectors.sam" 
        samtools_no_vector_orphans_sam_to_bam += " > " + vector_removal_folder + "orphans_no_vectors.bam"
        
        samtools_no_vector_orphans_bam_to_fastq = ">&2 echo samtools vr orphans pt 2 | " 
        samtools_no_vector_orphans_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -f 4" + " -0 "
        samtools_no_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_no_vectors.fastq "  
        samtools_no_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_no_vectors.bam"    
        
        samtools_vector_orphans_bam_to_fastq = ">&2 echo samtools vr orphans pt 3 | "
        samtools_vector_orphans_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -F 4" + " -0 " 
        samtools_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_vectors_only.fastq "  
        samtools_vector_orphans_bam_to_fastq += vector_removal_folder + "orphans_no_vectors.bam"
        
        bwa_vr_pair = ">&2 echo bwa vr pair | " 
        bwa_vr_pair += mpp.BWA + " mem -t " + self.Threads_str + " " 
        bwa_vr_pair += self.Vector_Contaminants + " " 
        bwa_vr_pair += blat_hr_folder + "pair_1_no_host.fastq " 
        bwa_vr_pair += blat_hr_folder + "pair_2_no_host.fastq" 
        bwa_vr_pair += " > " + vector_removal_folder + "pair_no_vectors.sam"
        
        samtools_vr_pair_sam_to_bam = ">&2 echo samtools vr pair pt 1 | "
        samtools_vr_pair_sam_to_bam += mpp.SAMTOOLS + " view -bS " 
        samtools_vr_pair_sam_to_bam += vector_removal_folder + "pair_no_vectors.sam" 
        samtools_vr_pair_sam_to_bam += " > " + vector_removal_folder + "pair_no_vectors.bam"
        
        samtools_no_vector_pair_bam_to_fastq = ">&2 echo samtools vr pair pt 2 | "
        samtools_no_vector_pair_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -f 13"
        samtools_no_vector_pair_bam_to_fastq += " -1 " + vector_removal_folder + "pair_1_no_vectors.fastq" 
        samtools_no_vector_pair_bam_to_fastq += " -2 " + vector_removal_folder + "pair_2_no_vectors.fastq " 
        samtools_no_vector_pair_bam_to_fastq += vector_removal_folder + "pair_no_vectors.bam"

        samtools_vector_pair_bam_to_fastq = ">&2 echo samtools vr pair pt 3 | "
        samtools_vector_pair_bam_to_fastq += mpp.SAMTOOLS + " fastq -n -F 4" 
        samtools_vector_pair_bam_to_fastq += " -1 " + vector_removal_folder + "pair_1_vectors_only.fastq" 
        samtools_vector_pair_bam_to_fastq += " -2 " + vector_removal_folder + "pair_2_vectors_only.fastq " 
        samtools_vector_pair_bam_to_fastq += vector_removal_folder + "pair_no_vectors.bam"
        
        make_blast_db_vector = ">&2 echo BLAST make db vectors | "
        make_blast_db_vector += mpp.Makeblastdb + " -in " + self.Vector_Contaminants + " -dbtype nucl"
        
        vsearch_filter_6 = ">&2 echo convert vr orphans for BLAT | " 
        vsearch_filter_6 += mpp.vsearch 
        vsearch_filter_6 += " --fastq_filter " + vector_removal_folder + "orphans_no_vectors.fastq" 
        vsearch_filter_6 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_6 += " --fastaout " + vector_removal_folder + "orphans_no_vectors.fasta"

        vsearch_filter_7 = ">&2 echo convert vr pair 1 for BLAT | "
        vsearch_filter_7 += mpp.vsearch 
        vsearch_filter_7 += " --fastq_filter " + vector_removal_folder + "pair_1_no_vectors.fastq" 
        vsearch_filter_7 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_7 += " --fastaout " + vector_removal_folder + "pair_1_no_vectors.fasta"
        
        vsearch_filter_8 = ">&2 echo convert vr pair 2 for BLAT | "
        vsearch_filter_8 += mpp.vsearch 
        vsearch_filter_8 += " --fastq_filter " + vector_removal_folder + "pair_2_no_vectors.fastq" 
        vsearch_filter_8 += " --fastq_ascii " + self.Qual_str 
        vsearch_filter_8 += " --fastaout " + vector_removal_folder + "pair_2_no_vectors.fasta"
        
        blat_vr_orphans = ">&2 echo BLAT vr orphans | "
        blat_vr_orphans += mpp.BLAT 
        blat_vr_orphans += " -noHead -minIdentity=90 -minScore=65 " 
        blat_vr_orphans += self.Vector_Contaminants + " " 
        blat_vr_orphans += vector_removal_folder + "orphans_no_vectors.fasta" 
        blat_vr_orphans += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " " 
        blat_vr_orphans += vector_removal_folder + "orphans_no_vectors.blatout"
        
        blat_vr_pair_1 = ">&2 echo BLAT vr pair 1 | "
        blat_vr_pair_1 += mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " 
        blat_vr_pair_1 += self.Vector_Contaminants + " " 
        blat_vr_pair_1 += vector_removal_folder + "pair_1_no_vectors.fasta" 
        blat_vr_pair_1 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " " 
        blat_vr_pair_1 += vector_removal_folder + "pair_1_no_vectors.blatout"
        
        blat_vr_pair_2 = ">&2 echo BLAT vr pair 2 | "
        blat_vr_pair_2 += mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " 
        blat_vr_pair_2 += self.Vector_Contaminants + " " 
        blat_vr_pair_2 += vector_removal_folder + "pair_2_no_vectors.fasta" 
        blat_vr_pair_2 += " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " " 
        blat_vr_pair_2 += vector_removal_folder + "pair_2_no_vectors.blatout"
        
        blat_containment_vector_folder = data_folder + "9_blat_containment_vr/"
        self.make_folder(blat_containment_vector_folder)
        blat_containment_vector_orphans = ">&2 echo BLAT contaminant orphans | " 
        blat_containment_vector_orphans += mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " 
        blat_containment_vector_orphans += vector_removal_folder + "orphans_no_vectors.fastq" + " " #in
        blat_containment_vector_orphans += vector_removal_folder + "orphans_no_vectors.blatout" + " " #in
        blat_containment_vector_orphans += blat_containment_vector_folder + "orphans_no_vectors.fastq" + " " #out
        blat_containment_vector_orphans += blat_containment_vector_folder + "orphans_vectors_only.fastq" #out
        
        blat_containment_vector_pair_1 = ">&2 echo BLAT contaminant pair 1 | "
        blat_containment_vector_pair_1 += mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " 
        blat_containment_vector_pair_1 += vector_removal_folder + "pair_1_no_vectors.fastq" + " " 
        blat_containment_vector_pair_1 += vector_removal_folder + "pair_1_no_vectors.blatout" + " " 
        blat_containment_vector_pair_1 += blat_containment_vector_folder + "pair_1_no_vectors.fastq" + " " 
        blat_containment_vector_pair_1 += blat_containment_vector_folder + "pair_1_vectors_only.fastq"
        
        blat_containment_vector_pair_2 = ">&2 echo BLAT contaminant pair 2 | "
        blat_containment_vector_pair_2 += mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " 
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
            # remove adapters
            adapter_removal_line,
            #trim things
            
            vsearch_merge,
            cat_glue,
            vsearch_filter_0,
            vsearch_filter_1,
            vsearch_filter_2,
            orphan_read_filter,
            
            cdhit_orphans,
            # # move_unpaired_cluster,
            cdhit_pair_1,
            cdhit_pair_2,
            # # move_paired_cluster,
            # #----host removal
            copy_host,
            bwa_hr_prep,
            # #----SAMTOOLS makes bam files
            samtools_hr_prep,
            bwa_hr_orphans,
            bwa_hr_pair,
            samtools_hr_orphans_sam_to_bam,
            samtools_no_host_orphans_bam_to_fastq,
            samtools_host_orphans_bam_to_fastq,
            samtools_host_pair_sam_to_bam,
            samtools_no_host_pair_bam_to_fastq,
            samtools_host_pair_bam_to_fastq,
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
            # # #-----vector removal
            copy_vector,
            bwa_vr_prep,
            samtools_vr_prep,
            bwa_vr_orphans,
            samtools_no_vector_orphans_sam_to_bam,
            samtools_no_vector_orphans_bam_to_fastq,
            samtools_vector_orphans_bam_to_fastq,
            bwa_vr_pair,
            samtools_vr_pair_sam_to_bam,
            samtools_no_vector_pair_bam_to_fastq,
            samtools_vector_pair_bam_to_fastq,
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
        #split, transform, and throw data into the rRNA filter.  Get back mRNA (goal) and rRNA (garbage)
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
        
        file_splitter_orphans = mpp.Python + " " + mpp.File_splitter + " " 
        file_splitter_orphans += dep_loc + "orphans_no_vectors.fastq " 
        file_splitter_orphans += orphan_split_folder + "orphans "
        file_splitter_orphans += str(file_split_count)
        
        file_splitter_pair_1 = mpp.Python + " " + mpp.File_splitter + " " 
        file_splitter_pair_1 += dep_loc + "pair_1_no_vectors.fastq " 
        file_splitter_pair_1 += pair_1_split_folder + "pair_1 "
        file_splitter_pair_1 += str(file_split_count)
        
        file_splitter_pair_2 = mpp.Python + " " + mpp.File_splitter + " " 
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
        
        convert_fastq_to_fasta = mpp.vsearch 
        convert_fastq_to_fasta += " --fastq_filter " + fastq_in
        convert_fastq_to_fasta += " --fastq_ascii " + self.Qual_str 
        convert_fastq_to_fasta += " --fastaout " + fasta_io
        #print("converting", fastq_in )
        #print("placing in", fasta_io )
        
        infernal_command = mpp.Infernal 
        infernal_command += " -o /dev/null --tblout " 
        infernal_command += infernal_out
        infernal_command += " --anytrunc --rfam -E 0.001 "
        infernal_command += mpp.Rfam + " "
        infernal_command += fasta_io
        
        rRNA_filtration = mpp.Python + " " 
        rRNA_filtration += mpp.rRNA_filter + " "
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
        COMMANDS_combine = []
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        final_folder = data_folder + "final_results/"
        preprocess_subfolder = os.getcwd() + "/" + preprocess_stage_name + "/"
        
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        self.make_folder(final_folder)
        
        
        if not(os.path.exists(preprocess_subfolder)):
            #start our own custom preprocess steps
            print("duplicate repopulation without preprocess not ready")
            
        else:
            #we ran a previous preprocess.  grab files
            #need 3, 5(clstr only), and mRNA from the 2nd stage.
            #for the mRNA, we don't really care if it is.  This stage is just supposed to add in the missing duplicates from something that was stripped.
            
            hq_path = preprocess_subfolder + "data/3_ar_quality_filter/"
            cluster_path = preprocess_subfolder + "data/5_remove_duplicates/"
            
            repop_orphans = ">&2 echo Duplication repopulate Orphans | "
            repop_orphans += mpp.Python + " " + mpp.duplicate_repopulate + " " 
            repop_orphans += hq_path + "orphans_hq.fastq" + " "   #in -> way back when things were quality-filtered.  
                                                                                            #      step 2 in preprocess.  could still contain rRNA
            repop_orphans += dep_loc + "mRNA/orphans.fastq" + " "      #in -> rRNA filtration output
            repop_orphans += cluster_path + "orphans_unique.fastq.clstr" + " "           #in -> duplicates filter output
            repop_orphans += final_folder + "orphans.fastq"        #out
            
            
            repop_pair_1 = ">&2 echo Duplication repopulate pair 1 | "
            repop_pair_1 += mpp.Python + " " + mpp.duplicate_repopulate + " " 
            repop_pair_1 += hq_path + "pair_1_hq.fastq" + " " 
            repop_pair_1 += dep_loc + "mRNA/pair_1.fastq" + " " 
            repop_pair_1 += cluster_path + "pair_1_unique.fastq.clstr" + " " 
            repop_pair_1 += final_folder + "pair_1.fastq"
            
            repop_pair_2 = ">&2 echo Duplication repopulate pair 2 | "
            repop_pair_2 += mpp.Python + " " + mpp.duplicate_repopulate + " " 
            repop_pair_2 += hq_path + "pair_2_hq.fastq" + " " 
            repop_pair_2 += dep_loc + "mRNA/pair_2.fastq" + " " 
            repop_pair_2 += cluster_path + "pair_2_unique.fastq.clstr" + " " 
            repop_pair_2 += final_folder + "pair_2.fastq"
        
            COMMANDS_Combine = [
            repop_orphans,
            repop_pair_1,
            repop_pair_2
            ]
        return COMMANDS_Combine

    def create_assemble_commands(self, stage_name, dependency_stage_name):
        subfolder = os.getcwd() + "/" + stage_name + "/"
        data_folder = subfolder + "data/"
        dep_loc = os.getcwd() + "/" + dependency_stage_name + "/data/final_results/"
        self.make_folder(subfolder)
        self.make_folder(data_folder)
        
        spades_folder = data_folder + "0_spades/"
        self.make_folder(spades_folder)
        
        #this assembles contigs
        spades = ">&2 echo Spades Contig assembly | "
        spades += mpp.Python + " " 
        spades += mpp.Spades + " -k 21,33,55,77 --meta" 
        spades += " -1 " + dep_loc + "pair_1.fastq" #in1 (pair 1)
        spades += " -2 " + dep_loc + "pair_2.fastq" #in2 (pair 2)
        spades += " -o " + spades_folder + "_SpadesOut" #out
        """
        bwa_index = mpp.BWA + " index -a bwtsw " + self.Contigs
        
        #calls BWA, then uses SAMTools to get a report
        bwa_pair_contigs = ">&2 echo BWA pair contigs | " 
        bwa_pair_contigs += mpp.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 " 
        bwa_pair_contigs += self.Contigs + " " 
        bwa_pair_contigs += self.Input_File1 + "_all_mRNA.fastq " 
        bwa_pair_contigs += self.Input_File2 + "_all_mRNA.fastq | " 
        bwa_pair_contigs += mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contig_paired.sam"
        
        bwa_orphans_contigs = ">&2 echo BWA orphan contigs | " 
        bwa_orphans_contigs += mpp.BWA + " mem -t " + self.Threads_str + " -B 40 -O 60 -E 10 -L 50 " 
        bwa_orphans_contigs += self.Contigs + " " 
        bwa_orphans_contigs += self.Input_Filepath + "_all_mRNA_unpaired.fastq | " 
        bwa_orphans_contigs += mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contig_unpaired.sam"
        
        contig_merge = ">&2 echo Contig merge | "
        contig_merge += mpp.Python + " " + mpp.Map_reads_contigs + " " 
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
                        #bwa_index,
                        #bwa_paired_contigs,
                        #bwa_unpaired_contigs#,
                        #contig_merge
                        ]
        return COMMANDS_Assemble
        
        
    def create_BWA_annotate_command(self):    
        
        bwa_contigs = mpp.BWA + " mem -t " + self.Threads_str + " " + mpp.DNA_DB + " " + self.Contigs + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contigs_BWA.sam"
        make_sam_1 = mpp.BWA + " mem -t " + self.Threads_str + " " + mpp.DNA_DB + " " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_unpaired_unmapped_BWA.sam"
        make_sam_2 = mpp.BWA + " mem -t " + self.Threads_str + " " + mpp.DNA_DB + " " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_paired_unmapped_BWA.sam"
        COMMANDS_Annotate_BWA = [
            bwa_contigs,
            make_sam_1,
            make_sam_2,
            mpp.Python + " " + mpp.Map_reads_gene_BWA + " " + mpp.DNA_DB + " " + self.Input_Filepath + "_contig_map.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + self.Contigs + " " + self.Input_Filepath + "_contigs_BWA.sam" + " " + self.Input_Filepath + "_contigs_n_BWA.fasta" + " " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " " + self.Input_Filepath + "_unpaired_unmapped_BWA.sam" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " " + self.Input_Filepath + "_paired_unmapped_BWA.sam" + " " + self.Input_File1 + "_unmapped_n_BWA.fasta" + " " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " " + self.Input_Filepath + "_paired_unmapped_BWA.sam" + " " + self.Input_File2 + "_unmapped_n_BWA.fasta"
        ]
        return COMMANDS_Annotate_BWA


        
    def create_BLAT_annotate_command(self, splits = 5):
        # example input filepath: self.Input_File1, self.Input_File2
        # example extension: [contigs_n_BWA, unpaired_unmapped_n_BWA, unmapped_BWA]-> leave out the .fasta.  it's implied
        # example datatype [contigs, unpaired, paired]
        self.Input_Filepath = input_filepath
        COMMANDS_Annotate_BLAT = []
        for i in range (1, splits+1):
            tag = "_" + str(i)
            blat_command = BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.DNA_DB_Prefix + tag + mpp.DNA_DB_Extension + " " + self.Input_Filepath + "_" + extension + ".fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads_str + " " + self.Input_Filepath + "_" + datatype + tag + ".blatout"
            COMMANDS_Annotate_BLAT.append(blat_command)
        
        final_cat = "cat " + self.Input_Filepath + "_" + datatype + "_[1-" + splits + "]" + ".blatout" + " > " + self.Input_Filepath + "_" + datatype + ".blatout"
        COMMANDS_Annotate_BLAT.append(final_cat)
        
        return COMMANDS_Annotate_BLAT

    def create_BLAT_pp_command(self):

        COMMANDS_Annotate_BLAT_Post = [
                        mpp.Python + " " + mpp.Map_reads_gene_BLAT + " " + mpp.DNA_DB + " " + self.Input_Filepath + "_contig_map.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + self.Input_Filepath + "_genes.fna" + " " + self.Input_Filepath + "_contigs_n_BWA.fasta" + " " + self.Input_Filepath + "_contigs" + ".blatout" + " " + self.Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + self.Input_Filepath + "_unpaired" + ".blatout" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File1 + "_unmapped_n_BWA.fasta" + " " + self.Input_File1 + "_paired" + ".blatout" + " " + self.Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File2 + "_unmapped_n_BWA.fasta" + " " + self.Input_File2 + "_paired" + ".blatout" + " " + self.Input_File2 + "_unmapped_n_BWA_BLAT.fasta"
                        ]
        return COMMANDS_Annotate_BLAT_Post

    def create_DIAMOND_annotate_command(self, count = 5):
        
        COMMANDS_Annotate_Diamond = []
        for i in range(1, count+1):
            tag = "_dmnd_tmp" + str(count)
            Diamond_command_list = ["mkdir -p " + os.path.splitext(self.Input_FName)[0] + tag,
                            mpp.DIAMOND + " blastx -p " + self.Threads_str + " -d " + mpp.Prot_DB + " -q " + self.Input_Filepath + "_contigs_n_BWA_BLAT" + ".fasta" + " -o " + self.Input_Filepath + "_contigs.dmdout" + " -f 6 -t " + os.path.splitext(self.Input_FName)[0] + tag + "-k 10 --id 85 --query-cover 65 --min-score 60"]
                            
            COMMANDS_Annotate_Diamond.append(Diamond_command_list)                
        return COMMANDS_Annotate_Diamond


        
    def create_DIAMOND_pp_command(self):    
        # the command just calls the merger program
        COMMANDS_Annotate_Diamond_Post = [
        mpp.Python + " " + mpp.Map_reads_prot_DMND + " " + mpp.Prot_DB + " " + self.Input_Filepath + "_contig_map.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + self.Input_Filepath + "_genes.fna" + " " + self.Input_Filepath + "_proteins.faa" + " " + self.Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + self.Input_Filepath + "_contigs.dmdout" + " " + self.Input_Filepath + "_contigs_n_BWA_BLAT_DMD.fasta" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_Filepath + "_unpaired.dmdout" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT_DMD.fasta" + " " + self.Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File1 + "_paired.dmdout" + " " + self.Input_File1 + "_unmapped_n_BWA_BLAT_DMD.fasta" + " " + self.Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File2 + "_paired.dmdout" + " " + self.Input_File2 + "_unmapped_n_BWA_BLAT_DMD.fasta"
        ]
        return COMMANDS_Annotate_Diamond_Post

    def create_taxonomic_annotation_command(self):
        
        
        get_taxa_from_gene = mpp.Python + " " + mpp.Annotated_taxid + " " + self.Input_Filepath + "_gene_map.tsv" + " " + mpp.accession2taxid + " " + self.Input_Filepath + "_TaxIDOut.tsv"
        kaiju_on_contigs = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Contigs + " -z " + self.Threads_str + " -o " + self.Input_Filepath + "_contigs_KaijuOut.tsv"
        
        kaiju_on_unpaired = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " -z " + self.Threads_str + " -o " + self.Input_Filepath + "_unpaired_KaijuOut.tsv"
        
        kaiju_on_paired = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " -j " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " -z " + self.Threads_str + " -o " + self.Input_Filepath + "_paired_KaijuOut.tsv"
        
        cat_kaiju = "cat " + self.Input_Filepath + "_contigs_KaijuOut.tsv" + " " + self.Input_Filepath + "_unpaired_KaijuOut.tsv" + " " + self.Input_Filepath + "_paired_KaijuOut.tsv" + " > " + self.Input_Filepath + "_KaijuOut.tsv"
        
        centrifuge_on_unmapped = mpp.Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -1 " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " -2 " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " -U " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + self.Qual_str + " -p " + self.Threads_str + " -S " + self.Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " --report-file " + self.Input_Filepath + "_unmapped_CentrifugeReport.txt"
        
        centrifuge_on_contigs = mpp.Centrifuge + " -f -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -U " + self.Contigs + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + self.Qual_str + " -p " + self.Threads_str + " -S " + self.Input_Filepath + "_contigs_CentrifugeOut.tsv" + " --report-file " + self.Input_Filepath + "_contigs_CentrifugeReport.txt"
        
        cat_centrifuge = "cat " + self.Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " " + self.Input_Filepath + "_contigs_CentrifugeOut.tsv" + " > " + self.Input_Filepath + "_CentrifugeOut.tsv"
        
        wevote_combine = mpp.Python + " " + mpp.Classification_combine + " " + self.Input_Filepath + "_contig_map.tsv" + " " + self.Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + self.Input_Filepath + "_TaxIDOut.tsv" + " " + self.Input_Filepath + "_KaijuOut.tsv" + " " + self.Input_Filepath + "_CentrifugeOut.tsv"
        
        wevote_call = mpp.WEVOTE + " -o " + self.Input_Filepath + "_WEVOTEOut" + " --db " + mpp.WEVOTEDB + " -c"
        
        awk_cleanup = "awk -F \'\\t\' \'{print \"C\\t\"$1\"\\t\"$9}\' " + os.path.join(self.Input_Filepath + "_WEVOTEOut", os.path.splitext(self.Input_FName)[0] + "_WEVOTEOut_WEVOTE_Details.txt") + " > " + self.Input_Filepath + "_WEVOTEOut.tsv"
        
        taxid_to_english = mpp.Python + " " + mpp.Contrain_classification + " " + "family" + " " + self.Input_Filepath + "_WEVOTEOut.tsv" + " " + mpp.Nodes + " " + mpp.Names + " " + self.Input_Filepath + "_WEVOTEOut_family.tsv"
        
        kaiju_to_krona = mpp.Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + self.Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt"
        
        awk_cleanup_krona = "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        
        kt_import_text_cleanup = mpp.ktImportText + " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        
        COMMANDS_Classify = [
            get_taxa_from_gene,
            kaiju_on_contigs,
            kaiju_on_unpaired,
            kaiju_on_paired,
            cat_kaiju,
            centrifuge_on_unmapped,
            centrifuge_on_contigs,
            cat_centrifuge,
            wevote_combine,
            "mkdir -p " + self.Input_Filepath + "_WEVOTEOut",
            "cp "  + self.Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + self.Input_Filepath + "_WEVOTEOut",
            "cd " + os.path.dirname(WEVOTE),
            wevote_call,
            "cd $PBS_O_WORKDIR",
            awk_cleanup,
            taxid_to_english,
            kaiju_to_krona,
            awk_cleanup_krona,
            kt_import_text_cleanup
            
            ]
        return COMMANDS_Classify 

    def create_EC_preprocess_command(self):    
        COMMANDS_EC_Preprocess = [
                        "mkdir -p " + self.EC_Split,
                        "mkdir -p " + self.EC_Output,
                        mpp.Python + " " + File_splitter + " " + "1000" + " " + self.Input_Filepath + "_proteins.faa" + " " + self.EC_Split,
                        ]
        return COMMANDS_EC_Preprocess

    def create_EC_detect_command(self, JobID_EC_preprocess):
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.EC_Split = os.path.join(self.Input_Filepath + "_EC_Annotation", "Split")
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        #This call has the qalter embedded.
        #get rid of this later
        self.Threads_str = Thread_count
        COMMANDS_Detect = [
                        "JOBS=$(" + mpp.Python + " " + mpp.Detect_Submit + " " + self.EC_Split + " " + self.EC_Output + " " + self.Threads_str + " " + JobID_EC_Preprocess.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
                        ]
        return COMMANDS_Detect

    def create_EC_detect_combine_command(self):
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.Input_FName = os.path.basename(Input_File)
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        COMMANDS_Combine_Detect = ["cat " + os.path.join(self.EC_Output, "Detect", "*.toppred") + " > " + os.path.join(self.EC_Output, "Detect", os.path.splitext(self.Input_FName)[0] + "_proteins.toppred")]

        return COMMANDS_Combine_Detect
        
    def create_EC_PRIAM_command(self): 
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        COMMANDS_PRIAM = [
        "mkdir -p " + os.path.join(self.EC_Output, "PRIAM"),
        "cd " + os.path.join(self.EC_Output, "PRIAM"),
        "java -jar" + " " + mpp.Priam + " -n " + os.path.splitext(self.Input_FName)[0] + "_PRIAM" + " -i " + self.Input_Filepath + "_proteins.faa" + " -p " + os.path.join(os.path.dirname(Priam), "PRIAM_MAR15") + " -od " + os.path.join(self.EC_Output, "PRIAM") +" -e T -pt 0.5 -mo -1 -mp 70 -cc T -cg T -bd " + mpp.BLAST_dir,
        ]
        return COMMANDS_PRIAM

    def create_EC_Diamond_command(self):    
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        COMMANDS_EC_Diamond = [
        "mkdir -p " + os.path.join(self.EC_Output, "Diamond"),
        "cd " + os.path.join(self.EC_Output, "Diamond"),
        mpp.DIAMOND + " blastp -p " + self.Threads_str + " --query "+ self.Input_Filepath + "_proteins.faa" + " --db "+ mpp.SWISS_PROT + " --outfmt "+ "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident" + " --out " + os.path.join(self.EC_Output, "Diamond", os.path.splitext(self.Input_FName)[0] + ".blastout") + " --evalue 0.0000000001 --max-target-seqs 1"
        ]
        return COMMANDS_EC_Diamond

    def create_EC_postprocess_command(self):
        COMMANDS_EC_Postprocess = [
        mpp.Python + " " + mpp.EC_Annotation_Post + " " + self.Input_Filepath + "_proteins.faa" + " " + self.EC_Output
        #Produce_Table
        ]
                
        return COMMANDS_EC_Postprocess

    def create_Network_generation_command(self):
        COMMANDS_Network = [
        mpp.Python + " " + mpp.RPKM + " " + mpp.Nodes + " " + self.Input_Filepath + "_WEVOTEOut.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + os.path.join(self.Input_Filepath + "_EC_Annotation", "Output", "Consolidated", os.path.splitext(self.Input_FName)[0] + "_proteins.ECs_All") + " " + self.Input_Filepath + "_RPKM.tsv" + " " + self.Input_Filepath + "_Cytoscape.tsv"
        ]
        return COMMANDS_Network

    def create_Join_command(self):
        wevote_cat = "cat " + self.Input_Filepath + "*/*_WEVOTEOut.tsv" + " > " + self.Input_Filepath + "_WEVOTEOut.tsv"
        taxid_to_english = mpp.Python + " " + Contrain_classification + " " + "family" + " " + self.Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + self.Input_Filepath + "_WEVOTEOut_family.tsv"
        kaiju_to_krona = Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + self.Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt"
        awk_cleanup_krona = "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + self.Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        kt_import_text_cleanup = ktImportText + " -o " + self.Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + self.Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        cat_gene_maps = "cat " + self.Input_Filepath + "*/*_gene_map.tsv" + " > " + self.Input_Filepath + "_gene_map.tsv"
        cat_EC = "cat " + os.path.join(self.Input_Filepath + "*/*_EC_Annotation", "Output", "Consolidated", "*_proteins.ECs_All") + " > " + self.Input_Filepath + "_proteins.ECs_All"
        get_RPKM = mpp.Python + " " + mpp.RPKM + " " + mpp.Nodes + " " + self.Input_Filepath + "_WEVOTEOut.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + self.Input_Filepath + "_proteins.ECs_All" + " " + self.Input_Filepath + "_RPKM.tsv" + " " + self.Input_Filepath + "_Cytoscape.tsv"
        
        COMMANDS_Join = [
        wevote_cat,
        taxid_to_english,
        kaiju_to_krona,
        awk_cleanup_krona,
        kt_import_text_cleanup,
        cat_gene_maps,
        cat_EC,
        get_RPKM
        ]
        
        return COMMANDS_Join
                        
                    
                    