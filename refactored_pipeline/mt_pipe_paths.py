import os

class tool_path_obj:
    def __init__ (self, mode = "scinet"):
            
        
        
        if(mode == "scinet"):
            #--------------------------------------------
            # reference paths
            script_path = "/home/j/jparkins/billyc59/parkinson_int_meta_pipe/"
            refactor_path = script_path + "refactored_pipeline/"
            UniVec_Core = "/home/j/jparkins/mobolaji/Databases/UniVec_Core.fasta"
            Adapter = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
        
            #----------------------------------------------------------
            # external tools
            
            self.cdhit_dup = "/home/j/jparkins/mobolaji/Tools/CDHIT/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup"
            self.Timmomatic = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
            self.AdapterRemoval = "/home/j/jparkins/mobolaji/Tools/AdapterRemoval/adapterremoval-2.2.0/build/AdapterRemoval"
            self.vsearch = "/home/j/jparkins/mobolaji/Tools/vsearch/vsearch-2.4.2-linux-x86_64/bin/vsearch"

            self.Flash = "/home/j/jparkins/mobolaji/Tools/Flash/FLASH-1.2.11/flash"
            self.BWA = "/home/j/jparkins/mobolaji/Tools/BWA/bwa-0.7.5a/bwa"
            self.SAMTOOLS = "/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/samtools"
            self.BLAT = "/home/j/jparkins/mobolaji/Tools/pBLAT/pblat/pblat"
            self.DIAMOND = "/home/j/jparkins/mobolaji/Tools/Diamond/diamond"
            self.Blastdbcmd = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/blastdbcmd"
            self.Makeblastdb = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/makeblastdb"

            self.Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmscan"
            self.Kaiju = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju"
            self.Kaiju2krona = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju2krona"
            self.ktImportText = "/home/j/jparkins/mobolaji/Tools/Krona/KronaTools-2.7/bin/ktImportText"
            self.Centrifuge = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge"
            self.Centrifuge_report = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge-kreport"
            self.kSLAM = "/home/j/jparkins/mobolaji/Tools/k-SLAM/k-SLAM/SLAM"
            self.Barrnap = "/home/j/jparkins/mobolaji/Tools/Barrnap/bin/barrnap"
            self.Priam = "/home/j/jparkins/mobolaji/Tools/PRIAM/PRIAM_search.jar"
            self.BLAST_dir = "/home/j/jparkins/mobolaji/Tools/BLAST/blast-2.2.26/bin/"
            self.WEVOTE = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/run_WEVOTE_PIPELINE.sh"
            self.WEVOTEDB = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/WEVOTEDB"



            
            #--------------------------------------------
            # custom scripts


            self.sam_trimmer = refactor_path + "sam_trimmer.py"
            self.sort_reads = refactor_path + "sort_reads_refactor.py"
            self.Perl = "/home/j/jparkins/mobolaji/perl"
            self.Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkins/mobolaji/python"
            self.DNA_DB = "/scratch/j/jparkins/mobolaji/Microbial_cds_db/microbial_all_cds.fasta"
            self.DNA_DB_Prefix = os.path.splitext(DNA_DB)[0]
            self.DNA_DB_Extension = os.path.splitext(DNA_DB)[1]
            self.Prot_DB = "/scratch/j/jparkins/mobolaji/NCBI_nr_db/nr"
            self.Host = "/home/j/jparkins/mobolaji/Databases/Mouse_cds.fasta"
            self.BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
            self.Fastqc = "/home/j/jparkins/mobolaji/Tools/FastQC/fastqc"
            self.Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"
            self.Filter_rRNA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate = refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs = script_path + "Map_read_contigs.py"
            self.orphaned_read_filter = refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter = refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter = refactor_path + "file_splitter.py"
            self.Sort_Reads = script_path + "Read_Classification/Sort_Reads.py"
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            self.rRNA_filter = refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA = script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT = script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND = script_path + "Map_read_prot_DMND.py"
            self.Spades = "/home/j/jparkins/mobolaji/Tools/SPAdes/SPAdes-3.9.1-Linux/bin/spades.py"

            self.EC_Annotation_Prep = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            self.Detect_Submit = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            self.EC_Annotation_Post = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            self.Detect = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"

            self.SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            self.Nodes = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
            self.Names = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
            self.Annotated_taxid = script_path + "Read_Classification/Get_TaxID.py"
            self.Contrain_classification = script_path + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine = script_path + "Read_Classification/Combine_WEVOTE.py"

            self.accession2taxid = "/scratch/j/jparkins/mobolaji/accession2taxid/accession2taxid"
            self.RPKM = script_path + "RPKM.py"
            
        elif(mode == "docker" or mode == "Docker"):
        
            script_path             = "/pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
            refactor_path           = script_path
            tool_path               = "/pipeline_tools/"
            
            #----------------------------------------------------------
            # Reference files
            # this is some NCBI reference file for Vectors
            self.UniVec_Core    = reference_file_path + "UniVec_Core.fasta"
            # this too.
            self.Adapter        = reference_file_path + "Trimmomatic_adapters/TruSeq3-PE-2.fa"
            self.Host           = reference_file_path + "Mouse_cds.fasta"
            self.Rfam           = reference_file_path + "Rfam.cm"
            self.DNA_DB         =   "/dump_site/microbial_all_cds.fasta"
            self.DNA_DB_Prefix  = os.path.splitext(self.DNA_DB)[0]
            self.DNA_DB_Extension = os.path.splitext(self.DNA_DB)[1]
            
            #----------------------------------------------------------
            # external tools
            
            self.cdhit_dup          = tool_path + "cdhit_dup/cd-hit-dup" 
            self.Timmomatic         = tool_path + "Trimmomatic/trimmomatic-0.36.jar"
            self.AdapterRemoval     = tool_path + "adapterremoval/AdapterRemoval"
            self.vsearch            = tool_path + "vsearch/vsearch"
            self.Flash              = tool_path + "FLASH/flash"
            self.BWA                = tool_path + "BWA/bwa"
            self.SAMTOOLS           = tool_path + "samtools/samtools"
            self.BLAT               = tool_path + "PBLAT/pblat"
            self.DIAMOND            = tool_path + "DIAMOND/diamond"
            self.Blastdbcmd         = tool_path + "BLAST_p/blastdbcmd"
            self.Makeblastdb        = tool_path + "BLAST_p/makeblastdb"
            self.Infernal           = tool_path + "infernal/cmscan"
            self.Kaiju              = tool_path + "kaiju/kaiju"
            self.Kaiju2krona        = tool_path + "kaiju/kaiju2krona"
            self.ktImportText       = tool_path + "KronaTools/scripts/ImportText.pl"
            self.Centrifuge         = tool_path + "centrifuge/centrifuge"
            self.Centrifuge_report  = tool_path + "centrifuge/centrifuge-kreport"
            self.kSLAM              = tool_path + "k-SLAM/SLAM"
            self.Barrnap            = tool_path + "barrnap/barrnap"
            self.Priam              = tool_path + "PRIAM_search/PRIAM_search.jar"
            self.BLAST_dir          = tool_path + "BLAST_p"
            self.WEVOTE             = tool_path + "WEVOTE/run_WEVOTE_PIPELINE.sh"
            self.WEVOTEDB           = tool_path + "WEVOTE/WEVOTEDB"
            self.Spades             = tool_path + "SPAdes/bin/spades.py"


            
            #--------------------------------------------
            # custom scripts

            self.map_read_contig_v2         = refactor_path + "map_read_contig_v2.py"
            self.sam_trimmer                = refactor_path + "sam_trimmer.py"
            self.contig_duplicate_remover   = script_path + "contig_duplicate_remover.py"
            self.sort_reads                 = refactor_path + "sort_reads_refactor.py"
            self.Filter_rRNA                = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate       = refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs          = script_path + "Map_read_contigs.py"
            self.orphaned_read_filter       = refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = refactor_path + "file_splitter.py"
            self.Sort_Reads                 = script_path + "Read_Classification/Sort_Reads.py"
            self.rRNA_filter                = refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = script_path + "Map_read_prot_DMND.py"
            
            self.EC_Annotation_Prep         = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            self.Detect_Submit              = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            self.EC_Annotation_Post         = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            self.Detect                     = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
            self.Annotated_taxid            = script_path + "Read_Classification/Get_TaxID.py"
            self.Constrain_classification    = script_path + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine     = script_path + "Read_Classification/Combine_WEVOTE.py"
            
            self.Perl = "/home/j/jparkins/mobolaji/perl"
            self.Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkins/mobolaji/python"
            
            self.Prot_DB = "/scratch/j/jparkins/mobolaji/NCBI_nr_db/nr"
            
            self.BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
            self.Fastqc = "/home/j/jparkins/mobolaji/Tools/FastQC/fastqc"
            
            
            
            
            
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            
            

            self.SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            self.Nodes = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
            self.Names = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
            

            self.accession2taxid = "/scratch/j/jparkins/mobolaji/accession2taxid/accession2taxid"
            self.RPKM = script_path + "RPKM.py"
        
        elif(mode == "singularity" or mode == "Singularity"):
            #temp place for scripts in refactored pipeline.  we'll move it once it's finished
            script_path             = "/home/j/jparkin/billyc59/parkinson_int_meta_pipe/refactored_pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
            refactor_path           = script_path
            tool_path               = "/pipeline_tools/"
            scratch_path            = "/scratch/j/jparkin/billyc59/"
            #----------------------------------------------------------
            # Reference files
            # this is some NCBI reference file for Vectors
            self.UniVec_Core    = reference_file_path + "UniVec_Core.fasta"
            # this too.
            self.Adapter        = reference_file_path + "Trimmomatic_adapters/TruSeq3-PE-2.fa"
            self.Host           = reference_file_path + "Mouse_cds.fasta"
            self.Rfam           = reference_file_path + "Rfam.cm"
            self.DNA_DB         = scratch_path +"microbial_all_cds.fasta"
            self.DNA_DB_Prefix  = os.path.splitext(self.DNA_DB)[0]
            self.DNA_DB_Extension = os.path.splitext(self.DNA_DB)[1]
            
            #----------------------------------------------------------
            # external tools
            
            self.cdhit_dup          = tool_path + "cdhit_dup/cd-hit-dup" 
            self.Timmomatic         = tool_path + "Trimmomatic/trimmomatic-0.36.jar"
            self.AdapterRemoval     = tool_path + "adapterremoval/AdapterRemoval"
            self.vsearch            = tool_path + "vsearch/vsearch"
            self.Flash              = tool_path + "FLASH/flash"
            self.BWA                = tool_path + "BWA/bwa"
            self.SAMTOOLS           = tool_path + "samtools/samtools"
            self.BLAT               = tool_path + "PBLAT/pblat"
            self.DIAMOND            = tool_path + "DIAMOND/diamond"
            self.Blastdbcmd         = tool_path + "BLAST_p/blastdbcmd"
            self.Makeblastdb        = tool_path + "BLAST_p/makeblastdb"
            self.Infernal           = tool_path + "infernal/cmscan"
            self.Kaiju              = tool_path + "kaiju/kaiju"
            self.Kaiju2krona        = tool_path + "kaiju/kaiju2krona"
            self.ktImportText       = tool_path + "KronaTools/scripts/ImportText.pl"
            self.Centrifuge         = tool_path + "centrifuge/centrifuge"
            self.Centrifuge_report  = tool_path + "centrifuge/centrifuge-kreport"
            self.kSLAM              = tool_path + "k-SLAM/SLAM"
            self.Barrnap            = tool_path + "barrnap/barrnap"
            self.Priam              = tool_path + "PRIAM_search/PRIAM_search.jar"
            self.BLAST_dir          = tool_path + "BLAST_p"
            self.WEVOTE             = tool_path + "WEVOTE/run_WEVOTE_PIPELINE.sh"
            self.WEVOTEDB           = tool_path + "WEVOTE/WEVOTEDB"
            self.Spades             = tool_path + "SPAdes/bin/spades.py"


            
            #--------------------------------------------
            # custom scripts

            self.map_read_contig_v2         = refactor_path + "map_read_contig_v2.py"
            self.sam_trimmer                = refactor_path + "sam_trimmer.py"
            self.contig_duplicate_remover   = script_path + "contig_duplicate_remover.py"
            self.sort_reads                 = refactor_path + "sort_reads_refactor.py"
            #self.Filter_rRNA                = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate       = refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs          = script_path + "Map_read_contigs.py"
            self.orphaned_read_filter       = refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = refactor_path + "file_splitter.py"
            #self.Sort_Reads                 = script_path + "Read_Classification/Sort_Reads.py"
            self.rRNA_filter                = refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = script_path + "Map_read_prot_DMND.py"
            
            self.RPKM = script_path + "RPKM.py"
            
            #self.EC_Annotation_Prep         = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            #self.Detect_Submit              = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            #self.EC_Annotation_Post         = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            #self.Detect                     = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
            #self.Annotated_taxid            = script_path + "Read_Classification/Get_TaxID.py"
            #self.Constrain_classification    = script_path + "Read_Classification/Constrain_Classification.py"
            #self.Classification_combine     = script_path + "Read_Classification/Combine_WEVOTE.py"
            
            #self.Perl = "/home/j/jparkins/mobolaji/perl"
            #self.Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkins/mobolaji/python"
            
            #self.Prot_DB = "/scratch/j/jparkins/mobolaji/NCBI_nr_db/nr"
            
            #self.BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
            #self.Fastqc = "/home/j/jparkins/mobolaji/Tools/FastQC/fastqc"
            
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            
            #self.SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            #self.Nodes = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
            #self.Names = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
            
            #self.accession2taxid = "/scratch/j/jparkins/mobolaji/accession2taxid/accession2taxid"
            
#from rRNA_Filter

#Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch"
#Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"

