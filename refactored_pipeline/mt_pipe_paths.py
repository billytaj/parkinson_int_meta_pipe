import os

class tool_path_obj:
    def __init__ (self, mode = "scinet"):
            
        
        
        if(mode == "scinet"):
            #--------------------------------------------
            # reference paths
            self.script_path = "/home/j/jparkin/billyc59/parkinson_int_meta_pipe/"
            self.refactor_path = self.script_path + "refactored_pipeline/"
            self.UniVec_Core = "/home/j/jparkin/mobolaji/Databases/UniVec_Core.fasta"
            self.Adapter = "/home/j/jparkin/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
        
            #----------------------------------------------------------
            # external tools
            
            self.cdhit_dup = "/home/j/jparkin/mobolaji/Tools/CDHIT/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup"
            self.Timmomatic = "/home/j/jparkin/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
            self.AdapterRemoval = "/home/j/jparkin/mobolaji/Tools/AdapterRemoval/adapterremoval-2.2.0/build/AdapterRemoval"
            self.vsearch = "/home/j/jparkin/mobolaji/Tools/vsearch/vsearch-2.4.2-linux-x86_64/bin/vsearch"

            self.Flash = "/home/j/jparkin/mobolaji/Tools/Flash/FLASH-1.2.11/flash"
            self.BWA = "/home/j/jparkin/mobolaji/Tools/BWA/bwa-0.7.5a/bwa"
            self.SAMTOOLS = "/home/j/jparkin/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/samtools"
            self.BLAT = "/home/j/jparkin/mobolaji/Tools/pBLAT/pblat/pblat"
            self.DIAMOND = "/home/j/jparkin/mobolaji/Tools/Diamond/diamond"
            self.Blastdbcmd = "/home/j/jparkin/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/blastdbcmd"
            self.Makeblastdb = "/home/j/jparkin/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/makeblastdb"

            self.Infernal = "/home/j/jparkin/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmscan"
            self.Kaiju = "/home/j/jparkin/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju"
            self.Kaiju2krona = "/home/j/jparkin/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju2krona"
            self.ktImportText = "/home/j/jparkin/mobolaji/Tools/Krona/KronaTools-2.7/bin/ktImportText"
            self.Centrifuge = "/home/j/jparkin/mobolaji/Tools/Centrifuge/centrifuge/centrifuge"
            self.Centrifuge_report = "/home/j/jparkin/mobolaji/Tools/Centrifuge/centrifuge/centrifuge-kreport"
            self.kSLAM = "/home/j/jparkin/mobolaji/Tools/k-SLAM/k-SLAM/SLAM"
            self.Barrnap = "/home/j/jparkin/mobolaji/Tools/Barrnap/bin/barrnap"
            self.Priam = "/home/j/jparkin/mobolaji/Tools/PRIAM/PRIAM_search.jar"
            self.BLAST_dir = "/home/j/jparkin/mobolaji/Tools/BLAST/blast-2.2.26/bin/"
            self.WEVOTE = "/home/j/jparkin/mobolaji/Tools/WEVOTE/WEVOTE/run_WEVOTE_PIPELINE.sh"
            self.WEVOTEDB = "/home/j/jparkin/mobolaji/Tools/WEVOTE/WEVOTE/WEVOTEDB"



            
            #--------------------------------------------
            # custom scripts


            self.sam_trimmer = self.refactor_path + "sam_trimmer.py"
            self.sort_reads = self.refactor_path + "sort_reads_refactor.py"
            self.Perl = "/home/j/jparkin/mobolaji/perl"
            self.Perl_Script_Dir = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkin/mobolaji/python"
            self.DNA_DB = "/scratch/j/jparkin/mobolaji/Microbial_cds_db/microbial_all_cds.fasta"
            self.DNA_DB_Prefix = os.path.splitext(self.DNA_DB)[0]
            self.DNA_DB_Extension = os.path.splitext(self.DNA_DB)[1]
            self.Prot_DB = "/scratch/j/jparkin/mobolaji/NCBI_nr_db/nr"
            self.Host = "/home/j/jparkin/mobolaji/Databases/Mouse_cds.fasta"
            self.BBMap_Dir = "/home/j/jparkin/mobolaji/Tools/BBMap/bbmap"
            self.Fastqc = "/home/j/jparkin/mobolaji/Tools/FastQC/fastqc"
            self.Rfam = "/home/j/jparkin/mobolaji/Databases/Rfam_rRNA.cm"
            self.Filter_rRNA = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate = self.refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs = self.script_path + "Map_read_contigs.py"
            self.orphaned_read_filter = self.refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter = self.refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter = self.refactor_path + "file_splitter.py"
            self.Sort_Reads = self.script_path + "Read_Classification/Sort_Reads.py"
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            self.rRNA_filter = self.refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA = self.script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT = self.script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND = self.script_path + "Map_read_prot_DMND.py"
            self.Spades = "/home/j/jparkin/mobolaji/Tools/SPAdes/SPAdes-3.9.1-Linux/bin/spades.py"

            self.EC_Annotation_Prep = self.script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            self.Detect_Submit = self.script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            self.EC_Annotation_Post = self.script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            self.Detect = "/home/j/jparkin/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"

            self.SWISS_PROT = "/home/j/jparkin/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            self.Nodes = "/home/j/jparkin/mobolaji/Databases/taxdump/nodes.dmp"
            self.Names = "/home/j/jparkin/mobolaji/Databases/taxdump/names.dmp"
            self.Annotated_taxid = self.script_path + "Read_Classification/Get_TaxID.py"
            self.Contrain_classification = self.script_path + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine = self.script_path + "Read_Classification/Combine_WEVOTE.py"

            self.accession2taxid = "/scratch/j/jparkin/mobolaji/accession2taxid/accession2taxid"
            self.RPKM = self.script_path + "RPKM.py"
            
        elif(mode == "docker" or mode == "Docker"):
        
            self.script_path             = "/pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
            self.refactor_path           = self.script_path
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

            self.map_read_contig_v2         = self.refactor_path + "map_read_contig_v2.py"
            self.sam_trimmer                = self.refactor_path + "sam_trimmer.py"
            self.contig_duplicate_remover   = self.script_path + "contig_duplicate_remover.py"
            self.sort_reads                 = self.refactor_path + "sort_reads_refactor.py"
            self.Filter_rRNA                = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate       = self.refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs          = self.script_path + "Map_read_contigs.py"
            self.orphaned_read_filter       = self.refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = self.refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = self.refactor_path + "file_splitter.py"
            self.Sort_Reads                 = self.script_path + "Read_Classification/Sort_Reads.py"
            self.rRNA_filter                = self.refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = self.script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = self.script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = self.script_path + "Map_read_prot_DMND.py"
            
            self.EC_Annotation_Prep         = self.script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            self.Detect_Submit              = self.script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            self.EC_Annotation_Post         = self.script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            self.Detect                     = "/home/j/jparkin/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
            self.Annotated_taxid            = self.script_path + "Read_Classification/Get_TaxID.py"
            self.Constrain_classification    = self.script_path + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine     = self.script_path + "Read_Classification/Combine_WEVOTE.py"
            
            self.Perl = "/home/j/jparkin/mobolaji/perl"
            self.Perl_Script_Dir = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkin/mobolaji/python"
            
            self.Prot_DB = "/scratch/j/jparkin/mobolaji/NCBI_nr_db/nr"
            
            self.BBMap_Dir = "/home/j/jparkin/mobolaji/Tools/BBMap/bbmap"
            self.Fastqc = "/home/j/jparkin/mobolaji/Tools/FastQC/fastqc"
            
            
            
            
            
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            
            

            self.SWISS_PROT = "/home/j/jparkin/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            self.Nodes = "/home/j/jparkin/mobolaji/Databases/taxdump/nodes.dmp"
            self.Names = "/home/j/jparkin/mobolaji/Databases/taxdump/names.dmp"
            

            self.accession2taxid = "/scratch/j/jparkin/mobolaji/accession2taxid/accession2taxid"
            self.RPKM = self.script_path + "RPKM.py"
        
        elif(mode == "singularity" or mode == "Singularity"):
            #temp place for scripts in refactored pipeline.  we'll move it once it's finished
            self.script_path             = "/home/j/jparkin//mobolaji/Metatranscriptome_Scripts/refactored_pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
            self.refactor_path           = self.script_path
            tool_path               = "/pipeline_tools/"
            scratch_path            = "/scratch/j/jparkin/billyc59/"
            project_path            = "/project/j/jparkin/Lab_Databases/"
            #----------------------------------------------------------
            # Reference files
            # this is some NCBI reference file for Vectors
            self.UniVec_Core        = reference_file_path + "UniVec_Core.fasta"
            # this too.
            self.Adapter            = reference_file_path + "Trimmomatic_adapters/TruSeq3-PE-2.fa"
            self.Host               = reference_file_path + "Mouse_cds.fasta"
            self.Rfam               = reference_file_path + "Rfam.cm"
            self.DNA_DB             = scratch_path +"Microbial_cds_db/microbial_all_cds.fasta"
            self.DNA_DB_Prefix      = os.path.splitext(self.DNA_DB)[0]
            self.DNA_DB_Extension   = os.path.splitext(self.DNA_DB)[1]
            self.Prot_DB            = project_path + "diamond_v0922/nr"
            
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

            self.map_read_contig_v2         = self.refactor_path + "map_read_contig_v2.py"
            self.sam_trimmer                = self.refactor_path + "sam_trimmer.py"
            self.contig_duplicate_remover   = self.script_path + "contig_duplicate_remover.py"
            self.sort_reads                 = self.refactor_path + "sort_reads_refactor.py"
            #self.Filter_rRNA                = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate       = self.refactor_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs          = self.script_path + "Map_read_contigs.py"
            self.orphaned_read_filter       = self.refactor_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = self.refactor_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = self.refactor_path + "file_splitter.py"
            #self.Sort_Reads                 = script_path + "Read_Classification/Sort_Reads.py"
            self.rRNA_filter                = self.refactor_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = self.script_path + "map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = self.script_path + "map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = self.script_path + "map_read_prot_DMND.py"
            
            self.RPKM = self.script_path + "RPKM.py"
            
            #self.EC_Annotation_Prep         = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            #self.Detect_Submit              = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            #self.EC_Annotation_Post         = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            #self.Detect                     = "/home/j/jparkin/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
            #self.Annotated_taxid            = script_path + "Read_Classification/Get_TaxID.py"
            #self.Constrain_classification    = script_path + "Read_Classification/Constrain_Classification.py"
            #self.Classification_combine     = script_path + "Read_Classification/Combine_WEVOTE.py"
            
            #self.Perl = "/home/j/jparkin/mobolaji/perl"
            #self.Perl_Script_Dir = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Xuejian"
            self.Python = "python3" #"/home/j/jparkin/mobolaji/python"
            
            #self.Prot_DB = "/scratch/j/jparkin/mobolaji/NCBI_nr_db/nr"
            
            #self.BBMap_Dir = "/home/j/jparkin/mobolaji/Tools/BBMap/bbmap"
            #self.Fastqc = "/home/j/jparkin/mobolaji/Tools/FastQC/fastqc"
            
            #rRNA_Split_Jobs = refactor_path + "rRNA_Split_Jobs.py"
            
            #self.SWISS_PROT = "/home/j/jparkin/mobolaji/Databases/uniprot_sprot_annotated.fasta"

            #self.Nodes = "/home/j/jparkin/mobolaji/Databases/taxdump/nodes.dmp"
            #self.Names = "/home/j/jparkin/mobolaji/Databases/taxdump/names.dmp"
            
            #self.accession2taxid = "/scratch/j/jparkin/mobolaji/accession2taxid/accession2taxid"
            
#from rRNA_Filter

#Infernal = "/home/j/jparkin/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmsearch"
#Rfam = "/home/j/jparkin/mobolaji/Databases/Rfam_rRNA.cm"

