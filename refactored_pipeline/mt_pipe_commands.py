import mt_pipe_paths as mpp


COMMANDS_Pre = [
                # remove adapters
                mpp.AdapterRemoval + " --file1 " + Input_File1 + ".fastq" + " --file2 " + Input_File2 + ".fastq" + " --qualitybase " + Qual + " --threads " + Threads + " --minlength " + "30" + " --basename " + os.path.splitext(Input_FName)[0] + "_AdapterRemoval" + " --trimqualities " + " --output1 " + Input_File1 + "_trimmed.fastq" + " --output2 " + Input_File2 + "_trimmed.fastq" + " --singleton " + Input_Filepath + "_singletons_trimmed.fastq",
                #trim things
                vsearch + " --fastq_mergepairs " + Input_File1 + "_trimmed.fastq" + " --reverse " + Input_File2 + "_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastqout " + Input_Filepath + "_overlap_trimmed.fastq" + " --fastqout_notmerged_fwd " + Input_File1 + "_paired_trimmed.fastq" + " --fastqout_notmerged_rev " + Input_File2 + "_paired_trimmed.fastq",
                "cat " + Input_Filepath + "_overlap_trimmed.fastq" + " " + Input_Filepath + "_singletons_trimmed.fastq" + " > " + Input_Filepath + "_unpaired_trimmed.fastq",
                vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_Filepath + "_unpaired_quality.fastq",
                vsearch + " --fastq_filter " + Input_File1 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File1 + "_quality.fastq",
                vsearch + " --fastq_filter " + Input_File2 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File2 + "_quality.fastq",
                Python + " " + Paired_Reads_Filter + " " + Input_File1 + "_quality.fastq" + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File2 + "_quality.fastq" + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_Filepath + "_unpaired_quality.fastq",
                cdhit_dup + " -i " + Input_Filepath + "_unpaired_quality.fastq" + " -o " + Input_Filepath + "_unpaired_unique.fastq",
                "mv " + Input_Filepath + "_unpaired_unique.fastq.clstr" + " " + Input_Filepath + "_unpaired.clstr",
                cdhit_dup + " -i " + Input_File1 + "_paired_quality.fastq" + " -i2 " + Input_File2 + "_paired_quality.fastq" + " -o " + Input_File1 + "_paired_unique.fastq" + " -o2 " + Input_File2 + "_paired_unique.fastq",
                "mv " + Input_File1 + "_paired_unique.fastq.clstr" + " " + Input_Filepath + "_paired.clstr",
                "cp " + Host + " " + Host_Contaminants,
                BWA + " index -a bwtsw " + Host_Contaminants,
                SAMTOOLS + " faidx " + Host_Contaminants,
                BWA + " mem -t " + Threads + " " + Host_Contaminants + " " + Input_Filepath + "_unpaired_unique.fastq" + " > " + Input_Filepath + "_unpaired_host_contaminants.sam",
                SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_host_contaminants.sam" + " > " + Input_Filepath + "_unpaired_host_contaminants.bam",
                SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam",
                SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam",
                BWA + " mem -t " + Threads + " " + Host_Contaminants + " " + Input_File1 + "_paired_unique.fastq" + " " + Input_File2 + "_paired_unique.fastq" + " > " + Input_Filepath + "_paired_host_contaminants.sam",
                SAMTOOLS + " view -bS " + Input_Filepath + "_paired_host_contaminants.sam" + " > " + Input_Filepath + "_paired_host_contaminants.bam",
                SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam",
                SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam",
                Makeblastdb + " -in " + Host_Contaminants + " -dbtype nucl",
                vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout",
                vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_host_contaminants.blatout",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_host_contaminants.fastq",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_host_contaminants.fastq",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_host_contaminants.fastq",
                "cp " + UniVec_Core + " " + Vector_Contaminants,
                BWA + " index -a bwtsw " + Vector_Contaminants,
                SAMTOOLS + " faidx " + Vector_Contaminants,
                BWA + " mem -t " + Threads + " " + Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_unpaired_vector_contaminants.sam",
                SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_vector_contaminants.sam" + " > " + Input_Filepath + "_unpaired_vector_contaminants.bam",
                SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam",
                SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam",
                BWA + " mem -t " + Threads + " " + Vector_Contaminants + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_paired_vector_contaminants.sam",
                SAMTOOLS + " view -bS " + Input_Filepath + "_paired_vector_contaminants.sam" + " > " + Input_Filepath + "_paired_vector_contaminants.bam",
                SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam",
                SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam",
                Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl",
                vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout",
                vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_vector_contaminants.blatout",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_vector_contaminants.fastq",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_vector_contaminants.fastq",
                Python + " " + BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_vector_contaminants.fastq",
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants",
                Python + " " + File_splitter + " " + "10000" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants",
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants",
                Python + " " + File_splitter + " " + "10000" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants",
                Python + " " + File_splitter + " " + "10000" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"
                ]
                
                
                
                
                
COMMANDS_rRNA = [
                "JOBS=$(" + Python + " " + rRNA_Split_Jobs + " " + Input_File + " " + JobID_Pre.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
                ]                
                

                COMMANDS_Combine = [
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_Filepath + "_mRNA_unpaired.fastq",
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_Filepath + "_rRNA_unpaired.fastq",
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File1 + "_mRNA.fastq",
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File1 + "_rRNA.fastq",
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File2 + "_mRNA.fastq",
                "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File2 + "_rRNA.fastq",
                Python + " " + Reduplicate + " " + Input_Filepath + "_unpaired_quality.fastq" + " " + Input_Filepath + "_mRNA_unpaired.fastq" + " " + Input_Filepath + "_unpaired.clstr" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq",
                Python + " " + Reduplicate + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File1 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File1 + "_all_mRNA.fastq",
                Python + " " + Reduplicate + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_File2 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File2 + "_all_mRNA.fastq"
                ]


COMMANDS_Assemble = [
                "mkdir -p " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0]) + "_SpadesOut",
                #Python + " " + Spades + " --rna" + " -1 " + Input_File1 + "_all_mRNA.fastq" + " -2 " + Input_File2 + "_all_mRNA.fastq" + " -s " + Input_Filepath + "_all_mRNA_unpaired.fastq" + " -o " + os.path.splitext(Input_FName)[0] + "_SpadesOut",
                Python + " " + Spades + " -k 21,33,55,77 --meta" + " -1 " + Input_File1 + "_all_mRNA.fastq" + " -2 " + Input_File2 + "_all_mRNA.fastq" + " -o " + os.path.splitext(Input_FName)[0] + "_SpadesOut",
                BWA + " index -a bwtsw " + Contigs,
                BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_File1 + "_all_mRNA.fastq " + Input_File2 + "_all_mRNA.fastq | " + SAMTOOLS + " view > " + Input_Filepath + "_contig_paired.sam",
                BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_Filepath + "_all_mRNA_unpaired.fastq | " + SAMTOOLS + " view > " + Input_Filepath + "_contig_unpaired.sam",
                Python + " " + Map_reads_contigs + " " + Input_File1 + "_all_mRNA.fastq" + " " + Input_File2 + "_all_mRNA.fastq" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq" + " " + Input_Filepath + "_contig_paired.sam" + " " + Input_Filepath + "_contig_unpaired.sam" + " " + Input_Filepath + "_contig_map.tsv"
                ]
COMMANDS_Annotate_BWA = [
                BWA + " mem -t " + Threads + " " + DNA_DB + " " + Contigs + " | " + SAMTOOLS + " view > " + Input_Filepath + "_contigs_BWA.sam",
                BWA + " mem -t " + Threads + " " + DNA_DB + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " | " + SAMTOOLS + " view > " + Input_Filepath + "_unpaired_unmapped_BWA.sam",
                BWA + " mem -t " + Threads + " " + DNA_DB + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " | " + SAMTOOLS + " view > " + Input_Filepath + "_paired_unmapped_BWA.sam",
                Python + " " + Map_reads_gene_BWA + " " + DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Contigs + " " + Input_Filepath + "_contigs_BWA.sam" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " " + Input_Filepath + "_unpaired_unmapped_BWA.sam" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File2 + "_unmapped_n_BWA.fasta",
                ]



COMMANDS_Annotate_BLAT1 = [
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_1" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_2" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_3" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_4" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_5" + ".blatout",
                "cat " + Input_Filepath + "_contigs" + "_[1-5]" + ".blatout" + " > " + Input_Filepath + "_contigs" + ".blatout"
                ]


COMMANDS_Annotate_BLAT2 = [
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_1" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_2" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_3" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_4" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_5" + ".blatout",
                "cat " + Input_Filepath + "_unpaired" + "_[1-5]" + ".blatout" + " > " + Input_Filepath + "_unpaired" + ".blatout"
                ]

COMMANDS_Annotate_BLAT3 = [
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_1" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_2" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_3" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_4" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_5" + ".blatout",
                "cat " + Input_File1 + "_paired" + "_[1-5]" + ".blatout" + " > " + Input_File1 + "_paired" + ".blatout"
                ]

COMMANDS_Annotate_BLAT4 = [
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_1" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_2" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_3" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_4" + ".blatout",
                BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_5" + ".blatout",
                "cat " + Input_File2 + "_paired" + "_[1-5]" + ".blatout" + " > " + Input_File2 + "_paired" + ".blatout"
                ]

COMMANDS_Annotate_BLAT_Post = [
                Python + " " + Map_reads_gene_BLAT + " " + DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_contigs" + ".blatout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_Filepath + "_unpaired" + ".blatout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File1 + "_paired" + ".blatout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_paired" + ".blatout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta"
                ]

COMMANDS_Annotate_Diamond1 = [
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp1",
                DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " -o " + Input_Filepath + "_contigs.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp1 -k 10 --id 85 --query-cover 65 --min-score 60",
                ]


COMMANDS_Annotate_Diamond2 = [
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp2",
                DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_Filepath + "_unpaired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp2 -k 10 --id 85 --query-cover 65 --min-score 60",
                ]

COMMANDS_Annotate_Diamond3 = [
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp3",
                DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_File1 + "_paired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp3 -k 10 --id 85 --query-cover 65 --min-score 60",
                ]

COMMANDS_Annotate_Diamond4 = [
                "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp4",
                DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_File2 + "_paired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp4 -k 10 --id 85 --query-cover 65 --min-score 60"
                ]


COMMANDS_Annotate_Diamond_Post = [
                Python + " " + Map_reads_prot_DMND + " " + Prot_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_proteins.faa" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_contigs.dmdout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT_DMD.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired.dmdout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_paired.dmdout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_paired.dmdout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT_DMD.fasta"
                ]


COMMANDS_Classify = [
                Python + " " + Annotated_taxid + " " + Input_Filepath + "_gene_map.tsv" + " " + accession2taxid + " " + Input_Filepath + "_TaxIDOut.tsv",
                Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Contigs + " -z " + Threads + " -o " + Input_Filepath + "_contigs_KaijuOut.tsv",
                Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_unpaired_KaijuOut.tsv",
                Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -j " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_paired_KaijuOut.tsv",
                "cat " + Input_Filepath + "_contigs_KaijuOut.tsv" + " " + Input_Filepath + "_unpaired_KaijuOut.tsv" + " " + Input_Filepath + "_paired_KaijuOut.tsv" + " > " + Input_Filepath + "_KaijuOut.tsv",
                Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -1 " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -2 " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -U " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_unmapped_CentrifugeReport.txt",
                Centrifuge + " -f -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -U " + Contigs + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_contigs_CentrifugeReport.txt",
                "cat " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " > " + Input_Filepath + "_CentrifugeOut.tsv",
                #kSLAM + " " + "--db=/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/" + " " + "--output-file=" + Input_Filepath + "_unpaired_kSLAMOut" + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq",
                #kSLAM + " " + "--db=/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/" + " " + "--output-file=" + Input_Filepath + "_paired_kSLAMOut" + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_File2 + "_all_mRNA_unmapped.fastq",
                #"sed \'s/^.*/C\\t&/\' " + Input_Filepath + "_unpaired_kSLAMOut_PerRead" + " " + Input_Filepath + "_paired_kSLAMOut_PerRead" + " > " + Input_Filepath + "_kSLAMOut.tsv",
                Python + " " + Classification_combine + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_TaxIDOut.tsv" + " " + Input_Filepath + "_KaijuOut.tsv" + " " + Input_Filepath + "_CentrifugeOut.tsv",# + " " + Input_Filepath + "_kSLAMOut.tsv",
                "mkdir -p " + Input_Filepath + "_WEVOTEOut",
                "cp "  + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_WEVOTEOut",
                "cd " + os.path.dirname(WEVOTE),
                WEVOTE + " -o " + Input_Filepath + "_WEVOTEOut" + " --db " + WEVOTEDB + " -c",
                "cd $PBS_O_WORKDIR",
                "awk -F \'\\t\' \'{print \"C\\t\"$1\"\\t\"$9}\' " + os.path.join(Input_Filepath + "_WEVOTEOut", os.path.splitext(Input_FName)[0] + "_WEVOTEOut_WEVOTE_Details.txt") + " > " + Input_Filepath + "_WEVOTEOut.tsv",
                Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
                Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
                "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
                ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
                ]
COMMANDS_EC_Preprocess = [
                "mkdir -p " + EC_Split,
                "mkdir -p " + EC_Output,
                Python + " " + File_splitter + " " + "1000" + " " + Input_Filepath + "_proteins.faa" + " " + EC_Split,
                ]


COMMANDS_Detect = [
                "JOBS=$(" + Python + " " + Detect_Submit + " " + EC_Split + " " + EC_Output + " " + Threads + " " + JobID_EC_Preprocess.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
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
                Python + " " + EC_Annotation_Post + " " + Input_Filepath + "_proteins.faa" + " " + EC_Output
                #Produce_Table
                ]


COMMANDS_Network = [
                Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + os.path.join(Input_Filepath + "_EC_Annotation", "Output", "Consolidated", os.path.splitext(Input_FName)[0] + "_proteins.ECs_All") + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
                ]


COMMANDS_Join = [
            "cat " + Input_Filepath + "*/*_WEVOTEOut.tsv" + " > " + Input_Filepath + "_WEVOTEOut.tsv",
            Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
            Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
            "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
            ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
            "cat " + Input_Filepath + "*/*_gene_map.tsv" + " > " + Input_Filepath + "_gene_map.tsv",
            "cat " + os.path.join(Input_Filepath + "*/*_EC_Annotation", "Output", "Consolidated", "*_proteins.ECs_All") + " > " + Input_Filepath + "_proteins.ECs_All",
            Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_proteins.ECs_All" + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
            ]
                
                