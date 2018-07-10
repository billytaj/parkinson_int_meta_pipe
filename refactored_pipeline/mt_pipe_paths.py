import os
import sys

class tool_path_obj:
    def __init__ (self, mode = "scinet", name = "None"):
            
        
        self.name = name
        
        if(mode == "docker" or mode == "Docker"):
        
            script_path             = "/pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
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

            self.map_read_contig_v2         = script_path + "map_read_contig_v2.py"
            self.sam_trimmer                = script_path + "sam_trimmer.py"
            self.contig_duplicate_remover   = script_path + "contig_duplicate_remover.py"
            self.sort_reads                 = script_path + "sort_reads_refactor.py"
            self.Filter_rRNA                = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
            self.duplicate_repopulate       = script_path +  "duplicate_repopulation.py"
            self.Map_reads_contigs          = script_path + "Map_read_contigs.py"
            self.orphaned_read_filter       = script_path + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = script_path + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = script_path + "seq_file_splitter.py"
            self.Sort_Reads                 = script_path + "Read_Classification/Sort_Reads.py"
            self.rRNA_filter                = script_path+"rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = script_path + "Map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = script_path + "Map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = script_path + "Map_read_prot_DMND.py"
            
            self.EC_Annotation_Prep         = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
            self.Detect_Submit              = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
            self.EC_Annotation_Post         = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
            self.Detect                     = "/home/j/jparkin/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
            self.Annotated_taxid            = script_path + "Read_Classification/Get_TaxID.py"
            self.Constrain_classification   = script_path + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine     = script_path + "Read_Classification/Combine_WEVOTE.py"
            
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
            if(self.name == "billy" or self.name == "Billy"):
                self.script_path        = "/home/j/jparkin/billyc59/parkinson_int_meta_pipe/refactored_pipeline/" 
            elif(name == "bj" or name == "BJ" or name == "Mobolaji" or name == "mobolaji"):
                self.script_path = "/home/j/jparkin/mobolaji/Metatranscriptome_Scripts/refactored_pipeline/"
            else:
                sys.exit("no name given.  unsure who to direct the tools to.  EXITING")
            #
            
            script_path             = self.script_path#"/pipeline/"
            reference_file_path     = "/pipeline_reference_files/"
            tool_path               = "/pipeline_tools/"
            scratch_path            = "/scratch/j/jparkin/billyc59/"
            database_path           = "/project/j/jparkin/Lab_Databases/"
            #----------------------------------------------------------
            # Reference files
            # this is some NCBI reference file for Vectors
            self.UniVec_Core        = database_path + "univec_core/UniVec_Core.fasta"
            # this too.
            self.Adapter            = database_path + "Trimmomatic_adapters/TruSeq3-PE-2.fa"
            self.Host               = database_path + "Mouse_cds/Mouse_cds.fasta"
            self.Rfam               = "/scratch/j/jparkin/billyc59/Rfam/Rfam.cm" #database_path + "Rfam.cm"
            self.DNA_DB             = database_path + "microbial_cds_db/microbial_all_cds.fasta"
            self.DNA_DB_Prefix      = os.path.splitext(self.DNA_DB)[0]
            self.DNA_DB_Extension   = os.path.splitext(self.DNA_DB)[1]
            self.Prot_DB            = database_path + "diamond_v0922/nr.dmnd"
            self.Prot_DB_plain      = database_path + "nr/nr"
            # currently not in singularity package
            self.accession2taxid    = "/scratch/j/jparkin/mobolaji/accession2taxid/accession2taxid"
            self.nodes              = database_path + "WEVOTE_db/nodes.dmp"
            self.names              = database_path + "WEVOTE_db/names.dmp"
            self.Kaiju_db           = database_path + "kaiju_db/kaiju_db_nr.fmi"
            self.Centrifuge_db      = database_path + "centrifuge_db/nt"
            self.SWISS_PROT         = database_path + "swiss_prot_db/swiss_prot_db"
            self.SWISS_PROT_map     = database_path + "swiss_prot_db/SwissProt_EC_Mapping.tsv"
            self.PriamDB            = database_path + "PRIAM_db/"
            self.DetectDB           = database_path + "DETECTv2"

            #----------------------------------------------------------
            # external tools
            self.Python             =             "python3"
            self.Python2            =             "python2"
            #Need to add
            self.Java               =             "java -jar"
            self.cdhit_dup          = tool_path + "cdhit_dup/cd-hit-dup" 
            self.Timmomatic         = tool_path + "Trimmomatic/trimmomatic-0.36.jar"
            self.AdapterRemoval     = tool_path + "adapterremoval/AdapterRemoval"
            self.vsearch            = tool_path + "vsearch/vsearch"
            self.Flash              = tool_path + "FLASH/flash"
            self.BWA                = tool_path + "BWA/bwa"
            self.SAMTOOLS           = tool_path + "samtools/samtools"
            self.BLAT               = tool_path + "PBLAT/pblat"
            self.DIAMOND            = tool_path + "DIAMOND/diamond"
            self.Blastp             = tool_path + "BLAST_p/blastp"
            self.Needle             = "/home/j/jparkin/mobolaji/Tools/EMBOSS/EMBOSS-6.6.0/emboss/needle"
            self.Blastdbcmd         = tool_path + "BLAST_p/blastdbcmd"
            self.Makeblastdb        = tool_path + "BLAST_p/makeblastdb"
            self.Infernal           = tool_path + "infernal/cmscan"
            self.Kaiju              = tool_path + "kaiju/kaiju"
            self.Kaiju2krona        = tool_path + "kaiju/kaiju2krona"
            self.ktImportText       = tool_path + "KronaTools/scripts/ImportText.pl"
            self.Centrifuge         = tool_path + "centrifuge/centrifuge"
            self.Centrifuge_report  = tool_path + "centrifuge/centrifuge-kreport"
            self.Centrifuge_build   = tool_path + "centrifuge-build"
            self.kSLAM              = tool_path + "k-SLAM/SLAM"
            self.Barrnap            = tool_path + "barrnap/barrnap"
            self.Priam              = tool_path + "PRIAM_search/PRIAM_search.jar"
            #Need to add
            self.Detect             = "/home/j/jparkin/mobolaji/Tools/DETECTv2/detect_2.01.py"#tool_path   + "detect_2_py3.py"
            self.BLAST_dir          = tool_path + "BLAST_p"
            self.WEVOTE             = tool_path + "WEVOTE/WEVOTE"
            self.WEVOTEDB           = tool_path + "WEVOTE/WEVOTE_db" #points to the location of taxdump, needed by WEVOTE
            self.Spades             = tool_path + "SPAdes/bin/spades.py"
            

            
            #--------------------------------------------
            # custom scripts

            self.map_read_contig_v2         = script_path   + "map_read_contig_v2.py"
            self.sam_trimmer                = script_path   + "sam_trimmer.py"
            self.contig_duplicate_remover   = script_path   + "contig_duplicate_remover.py"
            self.sort_reads                 = script_path   + "sort_reads_refactor.py"
            self.duplicate_repopulate       = script_path   + "duplicate_repopulation.py"
            self.orphaned_read_filter       = script_path   + "orphaned_pair_filter.py"
            self.BLAT_Contaminant_Filter    = script_path   + "BLAT_Contaminant_Filter.py"
            self.File_splitter              = script_path   + "seq_file_splitter.py"
            self.rRNA_filter                = script_path   + "rRNA_filter_v2.py"
            self.Map_reads_gene_BWA         = script_path   + "map_read_gene_BWA.py"
            self.Map_reads_gene_BLAT        = script_path   + "map_read_gene_BLAT.py"
            self.Map_reads_prot_DMND        = script_path   + "map_read_prot_DMND.py"
            
            self.RPKM                       = script_path   + "RPKM.py"

            #Will have to modify this script to take more explicit arguments
            self.EC_Annotation_Post         = script_path   + "EC_Consolidation.py"
            
            self.Annotated_taxid            = script_path   + "Read_Classification/Get_TaxID.py"
            self.Constrain_classification   = script_path   + "Read_Classification/Constrain_Classification.py"
            self.Classification_combine     = script_path   + "Read_Classification/Combine_WEVOTE.py"

            
            #self.Fastqc = "/home/j/jparkin/mobolaji/Tools/FastQC/fastqc"

