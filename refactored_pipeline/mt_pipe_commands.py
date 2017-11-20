#The functions here generate the pipeline commands.
#Each command module is made up of (often) many sub stages that are used to get the final result.
#If you want to move around the ordering, you'd do that here.

import os
import sys
import mt_pipe_paths as mpp

class mt_pipe_commands:
    def __init__ (self, Input_File, self.Quality_score, Thread_count):
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.Input_File1 = self.Input_Filepath + "1"
        self.Input_File2 = self.Input_Filepath + "2"
        self.Input_FName = os.path.basename(Input_File)
        
        self.Qual = self.Quality_score
        self.Input_Path = os.path.dirname(Input_File)
        self.self.Contigs = os.path.join(self.Input_Path, os.path.splitext(self.self.Input_FName)[0] + "_SpadesOut", "contigs.fasta")
        self.Threads = Thread_count
        self.EC_Split = os.path.join(self.Input_Filepath + "_EC_Annotation", "Split")
        self.EC_Output = os.path.join(self.Input_Filepath + "_EC_Annotation", "Output")
        

    def create_pre_command(self):

        Adapter_removal_line = mpp.AdapterRemoval + " --file1 " + self.Input_File1 + ".fastq" + " --file2 " + self.Input_File2 + ".fastq" + " --qualitybase " + self.Qual + " --threads " + self.Threads + " --minlength " + "30" + " --basename " + os.path.splitext(self.Input_FName)[0] + "_AdapterRemoval" + " --trimqualities " + " --output1 " + self.Input_File1 + "_trimmed.fastq" + " --output2 " + self.Input_File2 + "_trimmed.fastq" + " --singleton " + self.Input_Filepath + "_singletons_trimmed.fastq"
        
        Vsearch_merge = mpp.vsearch + " --fastq_mergepairs " + self.Input_File1 + "_trimmed.fastq" + " --reverse " + self.Input_File2 + "_trimmed.fastq" + " --fastq_ascii " + self.Qual + " --fastqout " + self.Input_Filepath + "_overlap_trimmed.fastq" + " --fastqout_notmerged_fwd " + self.Input_File1 + "_paired_trimmed.fastq" + " --fastqout_notmerged_rev " + self.Input_File2 + "_paired_trimmed.fastq"
        
        Cat_glue = "cat " + self.Input_Filepath + "_overlap_trimmed.fastq" + " " + self.Input_Filepath + "_singletons_trimmed.fastq" + " > " + self.Input_Filepath + "_unpaired_trimmed.fastq"
        vesearch_filter_0 = mpp.vsearch + " --fastq_filter " + self.Input_Filepath + "_unpaired_trimmed.fastq" + " --fastq_ascii " + self.Qual + " --fastq_maxee " + "2.0" + " --fastqout " + self.Input_Filepath + "_unpaired_quality.fastq"

        Vsearch_filter_1 = mpp.vsearch + " --fastq_filter " + self.Input_File1 + "_paired_trimmed.fastq" + " --fastq_ascii " + self.Qual + " --fastq_maxee " + "2.0" + " --fastqout " + self.Input_File1 + "_quality.fastq"
        
        Vsearch_filter_2 = mpp.vsearch + " --fastq_filter " + self.Input_File2 + "_paired_trimmed.fastq" + " --fastq_ascii " + self.Qual + " --fastq_maxee " + "2.0" + " --fastqout " + self.Input_File2 + "_quality.fastq"
        
        Paired_read_filter = mpp.Python + " " + mpp.Paired_Reads_Filter + " " + self.Input_File1 + "_quality.fastq" + " " + self.Input_File1 + "_paired_quality.fastq" + " " + self.Input_File2 + "_quality.fastq" + " " + self.Input_File2 + "_paired_quality.fastq" + " " + self.Input_Filepath + "_unpaired_quality.fastq"
        
        Cdhit_unpaired = mpp.cdhit_dup + " -i " + self.Input_Filepath + "_unpaired_quality.fastq" + " -o " + self.Input_Filepath + "_unpaired_unique.fastq"
        move_unpaired_cluster = "mv " + self.Input_Filepath + "_unpaired_unique.fastq.clstr" + " " + self.Input_Filepath + "_unpaired.clstr"
        
        Cdhit_paired = mpp.cdhit_dup + " -i " + self.Input_File1 + "_paired_quality.fastq" + " -i2 " + self.Input_File2 + "_paired_quality.fastq" + " -o " + self.Input_File1 + "_paired_unique.fastq" + " -o2 " + self.Input_File2 + "_paired_unique.fastq"
        
        move_paired_cluster = "mv " + self.Input_File1 + "_paired_unique.fastq.clstr" + " " + self.Input_Filepath + "_paired.clstr"
        
        copy_host = "cp " + mpp.Host + " " + mpp.Host_Contaminants
        
        bwa_host_remove_prep = mpp.BWA + " index -a bwtsw " + mpp.Host_Contaminants
        samtools_host_remove_prep = mpp.SAMTOOLS + " faidx " + mpp.Host_Contaminants
        bwa_host_remove = mpp.BWA + " mem -t " + self.Threads + " " + mpp.Host_Contaminants + " " + self.Input_Filepath + "_unpaired_unique.fastq" + " > " + self.Input_Filepath + "_unpaired_host_contaminants.sam"
        samtools_host_unpaired_sam_to_bam = mpp.SAMTOOLS + " view -bS " + self.Input_Filepath + "_unpaired_host_contaminants.sam" + " > " + self.Input_Filepath + "_unpaired_host_contaminants.bam"
        samtools_host_unmatched_unpaired_fastq_to_bam = mpp.SAMTOOLS + " fastq -n -f 4" + " -0 " + self.Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_host_contaminants.bam",
        samtools_host_unpaired_fastq_to_bam = mpp.SAMTOOLS + " fastq -n -F 4" + " -0 " + self.Input_Filepath + "_unpaired_host_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_host_contaminants.bam"                
        
        bwa_host_unpaired = mpp.BWA + " mem -t " + self.Threads + " " + mpp.Host_Contaminants + " " + self.Input_File1 + "_paired_unique.fastq" + " " + self.Input_File2 + "_paired_unique.fastq" + " > " + self.Input_Filepath + "_paired_host_contaminants.sam"
        
        samtools_host_paired_0 = mpp.SAMTOOLS + " view -bS " + self.Input_Filepath + "_paired_host_contaminants.sam" + " > " + self.Input_Filepath + "_paired_host_contaminants.bam"
        samtools_host_paired_1 = mpp.SAMTOOLS + " fastq -n -f 13" + " -1 " + self.Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " -2 " + self.Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + self.Input_Filepath + "_paired_host_contaminants.bam"
        samtools_host_paired_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -1 " + self.Input_File1 + "_paired_host_contaminants.fastq" + " -2 " + self.Input_File2 + "_paired_host_contaminants.fastq" + " " + self.Input_Filepath + "_paired_host_contaminants.bam"
        
        make_blast_db_host = mpp.Makeblastdb + " -in " + mpp.Host_Contaminants + " -dbtype nucl"
        
        vsearch_unmatched_unpaired = mpp.vsearch + " --fastq_filter " + self.Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + self.Qual + " --fastaout " + self.Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta"
        
        blat_unpaired_host_remove = mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Host_Contaminants + " " + self.Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads + " " + self.Input_Filepath + "_unpaired_n_host_contaminants.blatout"
        
        vsearch_unmatched_paired = mpp.vsearch + " --fastq_filter " + self.Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + self.Qual + " --fastaout " + self.Input_File1 + "_paired_n_BWA_host_contaminants.fasta"
        
        blat_paired_host_remove = mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Host_Contaminants + " " + self.Input_File1 + "_paired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads + " " + self.Input_File1 + "_paired_n_host_contaminants.blatout"
        
        blat_containment_host_unpaired = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_n_host_contaminants.blatout" + " " + self.Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_BLAT_host_contaminants.fastq"
        
        blat_containment_host_paired_1 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " " + self.Input_File1 + "_paired_n_host_contaminants.blatout" + " " + self.Input_File1 + "_paired_n_host_contaminants.fastq" + " " + self.Input_File1 + "_paired_BLAT_host_contaminants.fastq"
        
        blat_containment_host_paired_2 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + self.Input_File1 + "_paired_n_host_contaminants.blatout" + " " + self.Input_File2 + "_paired_n_host_contaminants.fastq" + " " + self.Input_File2 + "_paired_BLAT_host_contaminants.fastq"
        
        copy_vector = "cp " + mpp.UniVec_Core + " " + mpp.Vector_Contaminants
        
        bwa_vector_remove_prep = mpp.BWA + " index -a bwtsw " + mpp.Vector_Contaminants
        
        samtools_vector_remove_prep = mpp.SAMTOOLS + " faidx " + mpp.Vector_Contaminants
        
        bwa_vector_remove_unpaired = mpp.BWA + " mem -t " + self.Threads + " " + mpp.Vector_Contaminants + " " + self.Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " > " + self.Input_Filepath + "_unpaired_vector_contaminants.sam"
        
        samtools_unpaired_vector_remove_0 = mpp.SAMTOOLS + " view -bS " + self.Input_Filepath + "_unpaired_vector_contaminants.sam" + " > " + self.Input_Filepath + "_unpaired_vector_contaminants.bam"
        
        samtool_unpaired_vector_remove_1 = mpp.SAMTOOLS + " fastq -n -f 4" + " -0 " + self.Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_vector_contaminants.bam"    
        
        samtools_unpaired_vector_remove_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -0 " + self.Input_Filepath + "_unpaired_vector_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_vector_contaminants.bam"
        
        bwa_vector_remove_paired = mpp.BWA + " mem -t " + self.Threads + " " + mpp.Vector_Contaminants + " " + self.Input_File1 + "_paired_n_host_contaminants.fastq" + " " + self.Input_File2 + "_paired_n_host_contaminants.fastq" + " > " + self.Input_Filepath + "_paired_vector_contaminants.sam"
        
        samtools_paired_vector_remove_0 = mpp.SAMTOOLS + " view -bS " + self.Input_Filepath + "_paired_vector_contaminants.sam" + " > " + self.Input_Filepath + "_paired_vector_contaminants.bam"
        
        samtools_paired_vector_remove_1 = mpp.SAMTOOLS + " fastq -n -f 13" + " -1 " + self.Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " -2 " + self.Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + self.Input_Filepath + "_paired_vector_contaminants.bam"
        
        samtools_paired_vector_remove_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -1 " + self.Input_File1 + "_paired_vector_contaminants.fastq" + " -2 " + self.Input_File2 + "_paired_vector_contaminants.fastq" + " " + self.Input_Filepath + "_paired_vector_contaminants.bam"
        
        make_blast_db_vector = mpp.Makeblastdb + " -in " + mpp.Vector_Contaminants + " -dbtype nucl"
        
        vsearch_unpaired_vector_remove = mpp.vsearch + " --fastq_filter " + self.Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + self.Qual + " --fastaout " + self.Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta"
        
        blat_unpaired_vector_remove = BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Vector_Contaminants + " " + self.Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads + " " + self.Input_Filepath + "_unpaired_n_vector_contaminants.blatout"
        
        vsearch_paired_vector_remove = mpp.vsearch + " --fastq_filter " + self.Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + self.Qual + " --fastaout " + self.Input_File1 + "_paired_n_BWA_vector_contaminants.fasta"
        
        blat_paired_vector_remove = BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Vector_Contaminants + " " + self.Input_File1 + "_paired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads + " " + self.Input_File1 + "_paired_n_vector_contaminants.blatout"
        
        blat_containment_vector_unpaired = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_n_vector_contaminants.blatout" + " " + self.Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + self.Input_Filepath + "_unpaired_BLAT_vector_contaminants.fastq"
        blat_containment_vector_paired_1 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " " + self.Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + self.Input_File1 + "_paired_n_contaminants.fastq" + " " + self.Input_File1 + "_paired_BLAT_vector_contaminants.fastq"
        blat_containment_vector_paired_2 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + self.Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + self.Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + self.Input_File2 + "_paired_n_contaminants.fastq" + " " + self.Input_File2 + "_paired_BLAT_vector_contaminants.fastq"

        file_splitter_0 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + self.Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + os.path.splitext(self.Input_FName)[0] + "_unpaired_n_contaminants"
        
        file_splitter_1 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + self.Input_File1 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants"
        
        file_splitter_2 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + self.Input_File2 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants"
        
        COMMANDS_Pre = [
            # remove adapters
            adapter_removal_line,
            #trim things
            vsearch_merge,
            cat_glue,
            vsearch_filter_0,
            vsearch_filter_1,
            vsearch_filter_2,
            paired_read_filter,
            cdhit_unpaired,
            move_unpaired_cluster,
            cdhit_paired,
            move_paired_cluster,
            copy_host,
            bwa_host_remove_prep,
            # SAMTOOLS makes bam files
            samtools_host_remove,
            bwa_host_remove_unpaired,
            samtools_host_unpaired_sam_to_bam,
            samtools_host_unmatched_unpaired_fastq_to_bam,
            samtools_host_unpaired_fastq_to_bam,
            bwa_host_unpaired,
            samtools_host_paired_0,
            samtools_host_paired_1,
            samtools_host_paired_2,
            make_blast_db_host,
            vsearch_unmatched_unpaired,
            blat_unpaired_host_remove,
            vsearch_unmatched_paired,
            blat_paired_host_remove,
            blat_containment_host_unpaired,
            blat_containment_host_paired_1,
            blat_containment_host_paired_2,
            copy_vector,
            bwa_vector_remove_prep,
            samtools_vector_remove_prep,
            bwa_vector_remove_unpaired,
            samtools_unpaired_vector_remove_0,
            samtools_unpaired_vector_remove_1,
            samtools_unpaired_vector_remove_2,
            bwa_vector_remove_paired,
            samtools_paired_vector_remove_0,
            samtools_paired_vector_remove_1,
            samtools_paired_vector_remove_2,
            
            make_blast_db_vector, 
            
            vsearch_unpaired_vector_remove,
            blat_unpaired_vector_remove,
            vsearch_paired_vector_remove,
            blat_paired_vector_remove,
            
            blat_containment_vector_unpaired,
            blat_containment_vector_paired_1,
            blat_containment_vector_paired_2,
            
            
            "mkdir -p " + os.path.splitext(self.Input_FName)[0] + "_unpaired_n_contaminants",
            "mkdir -p " + os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants",
            file_splitter_0,
            file_splitter_1,
            file_splitter_2
        ]
        return COMMANDS_Pre            
                    
                    
                    
    def create_infernal_command(Input_File, qsub_job_id):                
        COMMANDS_rRNA = [
                        "JOBS=$(" + mpp.Python + " " + rRNA_Split_Jobs + " " + Input_File + " " + qsub_job_id.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
        ]                
        return COMMANDS_rRNA                

    def create_combine_command(Input_File):
        self.Input_Filepath = os.path.splitext(Input_File)[0]
        self.Input_File1 = self.Input_Filepath + "1"
        self.Input_File2 = self.Input_Filepath + "2"
        self.Input_FName = os.path.basename(Input_File)
        
        combine_unpaired_mrna = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(self.Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + self.Input_Filepath + "_mRNA_unpaired.fastq"

        combine_unpaired_rrna = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(self.Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + self.Input_Filepath + "_rRNA_unpaired.fastq"
        
        combine_pair_1_mrna_fastq = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants", os.path.basename(self.Input_File1) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + self.Input_File1 + "_mRNA.fastq"
        
        combine_pair_1_rrna_fastq = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants", os.path.basename(self.Input_File1) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + self.Input_File1 + "_rRNA.fastq"
        
        combine_pair_2_mrna_fastq = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants", os.path.basename(self.Input_File2) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + self.Input_File2 + "_mRNA.fastq"
        
        combine_pair_2_rrna_fastq = "cat " + os.path.join(os.path.splitext(self.Input_FName)[0] + "_paired_n_contaminants", os.path.basename(self.Input_File2) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + self.Input_File2 + "_rRNA.fastq"
        
        reduplicate_unpaired = mpp.Python + " " + mpp.Reduplicate + " " + self.Input_Filepath + "_unpaired_quality.fastq" + " " + self.Input_Filepath + "_mRNA_unpaired.fastq" + " " + self.Input_Filepath + "_unpaired.clstr" + " " + self.Input_Filepath + "_all_mRNA_unpaired.fastq"
        
        reduplicate_pair_1 = mpp.Python + " " + mpp.Reduplicate + " " + self.Input_File1 + "_paired_quality.fastq" + " " + self.Input_File1 + "_mRNA.fastq" + " " + self.Input_Filepath + "_paired.clstr" + " " + self.Input_File1 + "_all_mRNA.fastq"
        
        reduplicate_pair_2 = mpp.Python + " " + mpp.Reduplicate + " " + self.Input_File2 + "_paired_quality.fastq" + " " + self.Input_File2 + "_mRNA.fastq" + " " + self.Input_Filepath + "_paired.clstr" + " " + self.Input_File2 + "_all_mRNA.fastq"
        
        COMMANDS_Combine = [
        combine_unpaired_mrna,
        combine_unpaired_rrna,
        combine_pair_1_mrna_fastq,
        combine_pair_1_rrna_fastq,
        combine_pair_2_mrna_fastq,
        combine_pair_2_rrna_fastq,
        reduplicate_unpaired,
        reduplicate_pair_1,
        reduplicate_pair_2
        ]
        return COMMANDS_Combine

    def create_assemble_commands(self):
        #this assembles contigs
        spades = mpp.Python + " " + mpp.Spades + " -k 21,33,55,77 --meta" + " -1 " + self.Input_File1 + "_all_mRNA.fastq" + " -2 " + self.Input_File2 + "_all_mRNA.fastq" + " -o " + os.path.splitext(self.Input_FName)[0] + "_SpadesOut"
        
        bwa_index = mpp.BWA + " index -a bwtsw " + self.Contigs
        
        bwa_paired_contigs = mpp.BWA + " mem -t " + self.Threads + " -B 40 -O 60 -E 10 -L 50 " + self.Contigs + " " + self.Input_File1 + "_all_mRNA.fastq " + self.Input_File2 + "_all_mRNA.fastq | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contig_paired.sam"
        
        bwa_unpaired_contigs = mpp.BWA + " mem -t " + self.Threads + " -B 40 -O 60 -E 10 -L 50 " + self.Contigs + " " + self.Input_Filepath + "_all_mRNA_unpaired.fastq | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contig_unpaired.sam"
        
        contig_merge = mpp.Python + " " + mpp.Map_reads_contigs + " " + self.Input_File1 + "_all_mRNA.fastq" + " " + self.Input_File2 + "_all_mRNA.fastq" + " " + self.Input_Filepath + "_all_mRNA_unpaired.fastq" + " " + self.Input_Filepath + "_contig_paired.sam" + " " + self.Input_Filepath + "_contig_unpaired.sam" + " " + self.Input_Filepath + "_contig_map.tsv"
        
        COMMANDS_Assemble = [
                        "mkdir -p " + os.path.join(self.Input_Path, os.path.splitext(self.Input_FName)[0]) + "_SpadesOut",
                        spades,
                        bwa_index,
                        bwa_paired_contigs,
                        bwa_unpaired_contigs,
                        contig_merge
                        ]
        return COMMANDS_Assemble
        
        
    def create_BWA_annotate_command(self):    
        
        bwa_contigs = mpp.BWA + " mem -t " + self.Threads + " " + mpp.DNA_DB + " " + self.Contigs + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_contigs_BWA.sam"
        make_sam_1 = mpp.BWA + " mem -t " + self.Threads + " " + mpp.DNA_DB + " " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_unpaired_unmapped_BWA.sam"
        make_sam_2 = mpp.BWA + " mem -t " + self.Threads + " " + mpp.DNA_DB + " " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + self.Input_Filepath + "_paired_unmapped_BWA.sam"
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
            blat_command = BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.DNA_DB_Prefix + tag + mpp.DNA_DB_Extension + " " + self.Input_Filepath + "_" + extension + ".fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + self.Threads + " " + self.Input_Filepath + "_" + datatype + tag + ".blatout"
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
            Diamond_command_string = 
                            "mkdir -p " + os.path.splitext(self.Input_FName)[0] + tag,
                            mpp.DIAMOND + " blastx -p " + self.Threads + " -d " + mpp.Prot_DB + " -q " + self.Input_Filepath + "_contigs_n_BWA_BLAT" + ".fasta" + " -o " + self.Input_Filepath + "_contigs.dmdout" + " -f 6 -t " + os.path.splitext(self.Input_FName)[0] + tag + "-k 10 --id 85 --query-cover 65 --min-score 60",
                            
            COMMANDS_Annotate_Diamond.append(Diamond_command_string)                
        return COMMANDS_Annotate_Diamond


        
    def create_DIAMOND_pp_command(self):    
        # the command just calls the merger program
        COMMANDS_Annotate_Diamond_Post = [
        mpp.Python + " " + mpp.Map_reads_prot_DMND + " " + mpp.Prot_DB + " " + self.Input_Filepath + "_contig_map.tsv" + " " + self.Input_Filepath + "_gene_map.tsv" + " " + self.Input_Filepath + "_genes.fna" + " " + self.Input_Filepath + "_proteins.faa" + " " + self.Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + self.Input_Filepath + "_contigs.dmdout" + " " + self.Input_Filepath + "_contigs_n_BWA_BLAT_DMD.fasta" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_Filepath + "_unpaired.dmdout" + " " + self.Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT_DMD.fasta" + " " + self.Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File1 + "_paired.dmdout" + " " + self.Input_File1 + "_unmapped_n_BWA_BLAT_DMD.fasta" + " " + self.Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " " + self.Input_File2 + "_paired.dmdout" + " " + self.Input_File2 + "_unmapped_n_BWA_BLAT_DMD.fasta"
        ]
        return COMMANDS_Annotate_Diamond_Post

    def create_EC_classify_command(self):
        
        
        get_taxa_from_gene = mpp.Python + " " + mpp.Annotated_taxid + " " + self.Input_Filepath + "_gene_map.tsv" + " " + mpp.accession2taxid + " " + self.Input_Filepath + "_TaxIDOut.tsv"
        kaiju_on_contigs = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Contigs + " -z " + self.Threads + " -o " + self.Input_Filepath + "_contigs_KaijuOut.tsv"
        
        kaiju_on_unpaired = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " -z " + self.Threads + " -o " + self.Input_Filepath + "_unpaired_KaijuOut.tsv"
        
        kaiju_on_paired = mpp.Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " -j " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " -z " + self.Threads + " -o " + self.Input_Filepath + "_paired_KaijuOut.tsv"
        
        cat_kaiju = "cat " + self.Input_Filepath + "_contigs_KaijuOut.tsv" + " " + self.Input_Filepath + "_unpaired_KaijuOut.tsv" + " " + self.Input_Filepath + "_paired_KaijuOut.tsv" + " > " + self.Input_Filepath + "_KaijuOut.tsv"
        
        centrifuge_on_unmapped = mpp.Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -1 " + self.Input_File1 + "_all_mRNA_unmapped.fastq" + " -2 " + self.Input_File2 + "_all_mRNA_unmapped.fastq" + " -U " + self.Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + self.Qual + " -p " + self.Threads + " -S " + self.Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " --report-file " + self.Input_Filepath + "_unmapped_CentrifugeReport.txt"
        
        centrifuge_on_contigs = mpp.Centrifuge + " -f -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -U " + self.Contigs + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + self.Qual + " -p " + self.Threads + " -S " + self.Input_Filepath + "_contigs_CentrifugeOut.tsv" + " --report-file " + self.Input_Filepath + "_contigs_CentrifugeReport.txt"
        
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

    def create_EC_preprocess_command(self)    
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
        self.Threads = Thread_count
        COMMANDS_Detect = [
                        "JOBS=$(" + mpp.Python + " " + mpp.Detect_Submit + " " + self.EC_Split + " " + self.EC_Output + " " + self.Threads + " " + JobID_EC_Preprocess.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
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
        mpp.DIAMOND + " blastp -p " + self.Threads + " --query "+ self.Input_Filepath + "_proteins.faa" + " --db "+ mpp.SWISS_PROT + " --outfmt "+ "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident" + " --out " + os.path.join(self.EC_Output, "Diamond", os.path.splitext(self.Input_FName)[0] + ".blastout") + " --evalue 0.0000000001 --max-target-seqs 1"
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
                        
                    