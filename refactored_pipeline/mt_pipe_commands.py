#The functions here generate the pipeline commands.
#Each command module is made up of (often) many sub stages that are used to get the final result.
#If you want to move around the ordering, you'd do that here.

import os
import sys
import mt_pipe_paths as mpp

def create_pre_command(Input_File, Quality_score, Thread_count):
    Qual = Quality_score
    Threads = Thread_count
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    Input_FName = os.path.basename(Input_File)

    Adapter_removal_line = mpp.AdapterRemoval + " --file1 " + Input_File1 + ".fastq" + " --file2 " + Input_File2 + ".fastq" + " --qualitybase " + Qual + " --threads " + Threads + " --minlength " + "30" + " --basename " + os.path.splitext(Input_FName)[0] + "_AdapterRemoval" + " --trimqualities " + " --output1 " + Input_File1 + "_trimmed.fastq" + " --output2 " + Input_File2 + "_trimmed.fastq" + " --singleton " + Input_Filepath + "_singletons_trimmed.fastq"
    
    Vsearch_merge = mpp.vsearch + " --fastq_mergepairs " + Input_File1 + "_trimmed.fastq" + " --reverse " + Input_File2 + "_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastqout " + Input_Filepath + "_overlap_trimmed.fastq" + " --fastqout_notmerged_fwd " + Input_File1 + "_paired_trimmed.fastq" + " --fastqout_notmerged_rev " + Input_File2 + "_paired_trimmed.fastq"
    
    Cat_glue = "cat " + Input_Filepath + "_overlap_trimmed.fastq" + " " + Input_Filepath + "_singletons_trimmed.fastq" + " > " + Input_Filepath + "_unpaired_trimmed.fastq"
    vesearch_filter_0 = mpp.vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_Filepath + "_unpaired_quality.fastq"

    Vsearch_filter_1 = mpp.vsearch + " --fastq_filter " + Input_File1 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File1 + "_quality.fastq"
    
    Vsearch_filter_2 = mpp.vsearch + " --fastq_filter " + Input_File2 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File2 + "_quality.fastq"
    
    Paired_read_filter = mpp.Python + " " + mpp.Paired_Reads_Filter + " " + Input_File1 + "_quality.fastq" + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File2 + "_quality.fastq" + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_Filepath + "_unpaired_quality.fastq"
    
    Cdhit_unpaired = mpp.cdhit_dup + " -i " + Input_Filepath + "_unpaired_quality.fastq" + " -o " + Input_Filepath + "_unpaired_unique.fastq"
    move_unpaired_cluster = "mv " + Input_Filepath + "_unpaired_unique.fastq.clstr" + " " + Input_Filepath + "_unpaired.clstr"
    
    Cdhit_paired = mpp.cdhit_dup + " -i " + Input_File1 + "_paired_quality.fastq" + " -i2 " + Input_File2 + "_paired_quality.fastq" + " -o " + Input_File1 + "_paired_unique.fastq" + " -o2 " + Input_File2 + "_paired_unique.fastq"
    
    move_paired_cluster = "mv " + Input_File1 + "_paired_unique.fastq.clstr" + " " + Input_Filepath + "_paired.clstr"
    
    copy_host = "cp " + mpp.Host + " " + mpp.Host_Contaminants
    
    bwa_host_remove_prep = mpp.BWA + " index -a bwtsw " + mpp.Host_Contaminants
    samtools_host_remove_prep = mpp.SAMTOOLS + " faidx " + mpp.Host_Contaminants
    bwa_host_remove = mpp.BWA + " mem -t " + Threads + " " + mpp.Host_Contaminants + " " + Input_Filepath + "_unpaired_unique.fastq" + " > " + Input_Filepath + "_unpaired_host_contaminants.sam"
    samtools_host_unpaired_sam_to_bam = mpp.SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_host_contaminants.sam" + " > " + Input_Filepath + "_unpaired_host_contaminants.bam"
    samtools_host_unmatched_unpaired_fastq_to_bam = mpp.SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam",
    samtools_host_unpaired_fastq_to_bam = mpp.SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam"                
    
    bwa_host_unpaired = mpp.BWA + " mem -t " + Threads + " " + mpp.Host_Contaminants + " " + Input_File1 + "_paired_unique.fastq" + " " + Input_File2 + "_paired_unique.fastq" + " > " + Input_Filepath + "_paired_host_contaminants.sam"
    
    samtools_host_paired_0 = mpp.SAMTOOLS + " view -bS " + Input_Filepath + "_paired_host_contaminants.sam" + " > " + Input_Filepath + "_paired_host_contaminants.bam"
    samtools_host_paired_1 = mpp.SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam"
    samtools_host_paired_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam"
    
    make_blast_db_host = mpp.Makeblastdb + " -in " + mpp.Host_Contaminants + " -dbtype nucl"
    
    vsearch_unmatched_unpaired = mpp.vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta"
    
    blat_unpaired_host_remove = mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Host_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout"
    
    vsearch_unmatched_paired = mpp.vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta"
    
    blat_paired_host_remove = mpp.BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.Host_Contaminants + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_host_contaminants.blatout"
    
    blat_containment_host_unpaired = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_host_contaminants.fastq"
    
    blat_containment_host_paired_1 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_host_contaminants.fastq"
    
    blat_containment_host_paired_2 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_host_contaminants.fastq"
    
    copy_vector = "cp " + mpp.UniVec_Core + " " + mpp.Vector_Contaminants
    
    bwa_vector_remove_prep = mpp.BWA + " index -a bwtsw " + mpp.Vector_Contaminants
    
    samtools_vector_remove_prep = mpp.SAMTOOLS + " faidx " + mpp.Vector_Contaminants
    
    bwa_vector_remove_unpaired = mpp.BWA + " mem -t " + Threads + " " + mpp.Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_unpaired_vector_contaminants.sam"
    
    samtools_unpaired_vector_remove_0 = mpp.SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_vector_contaminants.sam" + " > " + Input_Filepath + "_unpaired_vector_contaminants.bam"
    
    samtool_unpaired_vector_remove_1 = mpp.SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam"    
    
    samtools_unpaired_vector_remove_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam"
    
    bwa_vector_remove_paired = mpp.BWA + " mem -t " + Threads + " " + mpp.Vector_Contaminants + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_paired_vector_contaminants.sam"
    
    samtools_paired_vector_remove_0 = mpp.SAMTOOLS + " view -bS " + Input_Filepath + "_paired_vector_contaminants.sam" + " > " + Input_Filepath + "_paired_vector_contaminants.bam"
    
    samtools_paired_vector_remove_1 = mpp.SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam"
    
    samtools_paired_vector_remove_2 = mpp.SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam"
    
    make_blast_db_vector = Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl"
    
    vsearch_unpaired_vector_remove = mpp.vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta"
    
    blat_unpaired_vector_remove = BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout"
    
    vsearch_paired_vector_remove = mpp.vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta"
    
    blat_paired_vector_remove = BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_vector_contaminants.blatout"
    
    blat_containment_vector_unpaired = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_vector_contaminants.fastq"
    blat_containment_vector_paired_1 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_vector_contaminants.fastq"
    blat_containment_vector_paired_2 = mpp.Python + " " + mpp.BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_vector_contaminants.fastq"

    file_splitter_0 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants"
    
    file_splitter_1 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"
    
    file_splitter_2 = mpp.Python + " " + mpp.File_splitter + " " + "10000" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"
    
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
        
        
        "mkdir -p " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants",
        "mkdir -p " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants",
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
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    Input_FName = os.path.basename(Input_File)
    
    combine_unpaired_mrna = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_Filepath + "_mRNA_unpaired.fastq"

    combine_unpaired_rrna = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_Filepath + "_rRNA_unpaired.fastq"
    
    combine_pair_1_mrna_fastq = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File1 + "_mRNA.fastq"
    
    combine_pair_1_rrna_fastq = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File1 + "_rRNA.fastq"
    
    combine_pair_2_mrna_fastq = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File2 + "_mRNA.fastq"
    
    combine_pair_2_rrna_fastq = "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File2 + "_rRNA.fastq"
    
    COMMANDS_Combine = [
    combine_unpaired_mrna,
    combine_unpaired_rrna,
    combine_pair_1_mrna_fastq,
    combine_pair_1_rrna_fastq,
    combine_pair_2_mrna_fastq,
    combine_pair_2_rrna_fastq,
    mpp.Python + " " + mpp.Reduplicate + " " + Input_Filepath + "_unpaired_quality.fastq" + " " + Input_Filepath + "_mRNA_unpaired.fastq" + " " + Input_Filepath + "_unpaired.clstr" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq",
    mpp.Python + " " + mpp.Reduplicate + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File1 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File1 + "_all_mRNA.fastq",
    mpp.Python + " " + mpp.Reduplicate + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_File2 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File2 + "_all_mRNA.fastq"
    ]
    return COMMANDS_Combine

def create_assemble_commands(Input_File, Thread_count, Contigs):
    #this assembles contigs
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    Input_Path = os.path.dirname(Input_File)
    Input_FName = os.path.basename(Input_File)
    Threads = Thread_count
    
    spades = mpp.Python + " " + mpp.Spades + " -k 21,33,55,77 --meta" + " -1 " + Input_File1 + "_all_mRNA.fastq" + " -2 " + Input_File2 + "_all_mRNA.fastq" + " -o " + os.path.splitext(Input_FName)[0] + "_SpadesOut"
    
    bwa_index = mpp.BWA + " index -a bwtsw " + Contigs
    
    bwa_paired_contigs = mpp.BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_File1 + "_all_mRNA.fastq " + Input_File2 + "_all_mRNA.fastq | " + mpp.SAMTOOLS + " view > " + Input_Filepath + "_contig_paired.sam"
    
    bwa_unpaired_contigs = mpp.BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_Filepath + "_all_mRNA_unpaired.fastq | " + mpp.SAMTOOLS + " view > " + Input_Filepath + "_contig_unpaired.sam"
    
    contig_merge = mpp.Python + " " + mpp.Map_reads_contigs + " " + Input_File1 + "_all_mRNA.fastq" + " " + Input_File2 + "_all_mRNA.fastq" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq" + " " + Input_Filepath + "_contig_paired.sam" + " " + Input_Filepath + "_contig_unpaired.sam" + " " + Input_Filepath + "_contig_map.tsv"
    
    COMMANDS_Assemble = [
                    "mkdir -p " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0]) + "_SpadesOut",
                    spades,
                    bwa_index,
                    bwa_paired_contigs,
                    bwa_unpaired_contigs,
                    contig_merge
                    ]
    return COMMANDS_Assemble
    
    
def create_BWA_annotate_command(Contigs):    
    
    bwa_contigs = mpp.BWA + " mem -t " + Threads + " " + mpp.DNA_DB + " " + Contigs + " | " + mpp.SAMTOOLS + " view > " + Input_Filepath + "_contigs_BWA.sam"
    make_sam_1 = mpp.BWA + " mem -t " + Threads + " " + mpp.DNA_DB + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + Input_Filepath + "_unpaired_unmapped_BWA.sam"
    make_sam_2 = mpp.BWA + " mem -t " + Threads + " " + mpp.DNA_DB + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " | " + mpp.SAMTOOLS + " view > " + Input_Filepath + "_paired_unmapped_BWA.sam"
    COMMANDS_Annotate_BWA = [
        bwa_contigs,
        make_sam_1,
        make_sam_2,
        mpp.Python + " " + mpp.Map_reads_gene_BWA + " " + mpp.DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Contigs + " " + Input_Filepath + "_contigs_BWA.sam" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " " + Input_Filepath + "_unpaired_unmapped_BWA.sam" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File2 + "_unmapped_n_BWA.fasta"
    ]
    return COMMANDS_Annotate_BWA


    
def create_BLAT_annotate_command(input_filepath, extension, datatype, splits = 5):
    # example input filepath: Input_File1, Input_File2
    # example extension: [contigs_n_BWA, unpaired_unmapped_n_BWA, unmapped_BWA]-> leave out the .fasta.  it's implied
    # example datatype [contigs, unpaired, paired]
    Input_Filepath = input_filepath
    COMMANDS_Annotate_BLAT = []
    for i in range (1, splits+1):
        tag = "_" + str(i)
        blat_command = BLAT + " -noHead -minIdentity=90 -minScore=65 " + mpp.DNA_DB_Prefix + tag + mpp.DNA_DB_Extension + " " + Input_Filepath + "_" + extension + ".fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_" + datatype + tag + ".blatout"
        COMMANDS_Annotate_BLAT.append(blat_command)
    
    final_cat = "cat " + Input_Filepath + "_" + datatype + "_[1-" + splits + "]" + ".blatout" + " > " + Input_Filepath + "_" + datatype + ".blatout"
    COMMANDS_Annotate_BLAT.append(final_cat)
    
    return COMMANDS_Annotate_BLAT

def create_BLAT_pp_command(Input_File):

    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    Input_Path = os.path.dirname(Input_File)
    COMMANDS_Annotate_BLAT_Post = [
                    mpp.Python + " " + mpp.Map_reads_gene_BLAT + " " + mpp.DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_contigs" + ".blatout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_Filepath + "_unpaired" + ".blatout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File1 + "_paired" + ".blatout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_paired" + ".blatout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta"
                    ]
    return COMMANDS_Annotate_BLAT_Post

def create_DIAMOND_annotate_command(Input_File, datatype, Thread_count, count = 5):
    Threads = Thread_count
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_FName = os.path.basename(Input_File)
    COMMANDS_Annotate_Diamond = []
    for i in range(1, count+1):
        tag = "_dmnd_tmp" + str(count)
        Diamond_command_string = 
                        "mkdir -p " + os.path.splitext(Input_FName)[0] + tag,
                        mpp.DIAMOND + " blastx -p " + Threads + " -d " + mpp.Prot_DB + " -q " + Input_Filepath + "_contigs_n_BWA_BLAT" + ".fasta" + " -o " + Input_Filepath + "_contigs.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + tag + "-k 10 --id 85 --query-cover 65 --min-score 60",
                        
        COMMANDS_Annotate_Diamond.append(Diamond_command_string)                
    return COMMANDS_Annotate_Diamond


    
def create_DIAMOND_pp_command(Input_File):    
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_FName = os.path.basename(Input_File)
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    COMMANDS_Annotate_Diamond_Post = [
    mpp.Python + " " + mpp.Map_reads_prot_DMND + " " + mpp.Prot_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_proteins.faa" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_contigs.dmdout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT_DMD.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired.dmdout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_paired.dmdout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_paired.dmdout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT_DMD.fasta"
    ]
    return COMMANDS_Annotate_Diamond_Post

def create_EC_classify_command(Input_File, Thread_count, Quality_score):
    Qual = Quality_score
    Input_Path = os.path.dirname(Input_File)
    Contigs = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_SpadesOut", "contigs.fasta")
    Threads = Thread_count
    Input_Filepath = os.path.splitext(Input_File)[0]
    Input_FName = os.path.basename(Input_File)
    Input_File1 = Input_Filepath + "1"
    Input_File2 = Input_Filepath + "2"
    
    get_taxa_from_gene = mpp.Python + " " + Annotated_taxid + " " + Input_Filepath + "_gene_map.tsv" + " " + mpp.accession2taxid + " " + Input_Filepath + "_TaxIDOut.tsv"
    kaiju_on_contigs = Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Contigs + " -z " + Threads + " -o " + Input_Filepath + "_contigs_KaijuOut.tsv"
    
    kaiju_on_unpaired = Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_unpaired_KaijuOut.tsv"
    
    kaiju_on_paired = Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -j " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_paired_KaijuOut.tsv"
    
    cat_kaiju = "cat " + Input_Filepath + "_contigs_KaijuOut.tsv" + " " + Input_Filepath + "_unpaired_KaijuOut.tsv" + " " + Input_Filepath + "_paired_KaijuOut.tsv" + " > " + Input_Filepath + "_KaijuOut.tsv"
    
    centrifuge_on_unmapped = mpp.Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -1 " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -2 " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -U " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_unmapped_CentrifugeReport.txt"
    
    centrifuge_on_contigs = mpp.Centrifuge + " -f -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -U " + Contigs + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_contigs_CentrifugeReport.txt"
    
    cat_centrifuge = "cat " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " > " + Input_Filepath + "_CentrifugeOut.tsv"
    
    COMMANDS_Classify = [
        get_taxa_from_gene,
        kaiju_on_contigs,
        kaiju_on_unpaired,
        kaiju_on_paired,
        cat_kaiju,
        centrifuge_on_unmapped,
        centrifuge_on_contigs,
        cat_centrifuge,
        
        mpp.Python + " " + Classification_combine + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_TaxIDOut.tsv" + " " + Input_Filepath + "_KaijuOut.tsv" + " " + Input_Filepath + "_CentrifugeOut.tsv",
        
        "mkdir -p " + Input_Filepath + "_WEVOTEOut",
        
        "cp "  + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_WEVOTEOut",
        
        "cd " + os.path.dirname(WEVOTE),
        
        WEVOTE + " -o " + Input_Filepath + "_WEVOTEOut" + " --db " + WEVOTEDB + " -c",
        
        "cd $PBS_O_WORKDIR",
        
        "awk -F \'\\t\' \'{print \"C\\t\"$1\"\\t\"$9}\' " + os.path.join(Input_Filepath + "_WEVOTEOut", os.path.splitext(Input_FName)[0] + "_WEVOTEOut_WEVOTE_Details.txt") + " > " + Input_Filepath + "_WEVOTEOut.tsv",
        
        mpp.Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
        
        Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
        
        "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
        
        ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
        
        ]
COMMANDS_EC_Preprocess = [
                "mkdir -p " + EC_Split,
                "mkdir -p " + EC_Output,
                mpp.Python + " " + File_splitter + " " + "1000" + " " + Input_Filepath + "_proteins.faa" + " " + EC_Split,
                ]


COMMANDS_Detect = [
                "JOBS=$(" + mpp.Python + " " + Detect_Submit + " " + EC_Split + " " + EC_Output + " " + Threads + " " + JobID_EC_Preprocess.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
                ]

COMMANDS_Combine_Detect = ["cat " + os.path.join(EC_Output, "Detect", "*.toppred") + " > " + os.path.join(EC_Output, "Detect", os.path.splitext(Input_FName)[0] + "_proteins.toppred")]


COMMANDS_PRIAM = [
                "mkdir -p " + os.path.join(EC_Output, "PRIAM"),
                "cd " + os.path.join(EC_Output, "PRIAM"),
                "java -jar" + " " + Priam + " -n " + os.path.splitext(Input_FName)[0] + "_PRIAM" + " -i " + Input_Filepath + "_proteins.faa" + " -p " + os.path.join(os.path.dirname(Priam), "PRIAM_MAR15") + " -od " + os.path.join(EC_Output, "PRIAM") +" -e T -pt 0.5 -mo -1 -mp 70 -cc T -cg T -bd " + BLAST_dir,
                ]


COMMANDS_EC_Diamond = [
                "mkdir -p " + os.path.join(EC_Output, "Diamond"),
                "cd " + os.path.join(EC_Output, "Diamond"),
                DIAMOND + " blastp -p " + Threads + " --query "+ Input_Filepath + "_proteins.faa" + " --db "+ SWISS_PROT + " --outfmt "+ "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident" + " --out " + os.path.join(EC_Output, "Diamond", os.path.splitext(Input_FName)[0] + ".blastout") + " --evalue 0.0000000001 --max-target-seqs 1"
                ]


COMMANDS_EC_Postprocess = [
                mpp.Python + " " + EC_Annotation_Post + " " + Input_Filepath + "_proteins.faa" + " " + EC_Output
                #Produce_Table
                ]


COMMANDS_Network = [
                mpp.Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + os.path.join(Input_Filepath + "_EC_Annotation", "Output", "Consolidated", os.path.splitext(Input_FName)[0] + "_proteins.ECs_All") + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
                ]


COMMANDS_Join = [
            "cat " + Input_Filepath + "*/*_WEVOTEOut.tsv" + " > " + Input_Filepath + "_WEVOTEOut.tsv",
            mpp.Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
            Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
            "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
            ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
            "cat " + Input_Filepath + "*/*_gene_map.tsv" + " > " + Input_Filepath + "_gene_map.tsv",
            "cat " + os.path.join(Input_Filepath + "*/*_EC_Annotation", "Output", "Consolidated", "*_proteins.ECs_All") + " > " + Input_Filepath + "_proteins.ECs_All",
            mpp.Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_proteins.ECs_All" + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
            ]
                
                