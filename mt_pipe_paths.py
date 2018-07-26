import os
import sys

class tool_path_obj:
    def __init__ (self, config):

        self.config = config

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
        self.Rfam               = database_path + "Rfam/Rfam.cm"
        self.DNA_DB             = database_path + "ChocoPhlAn/ChocoPhlAn.fasta"#"microbial_cds_db/microbial_all_cds.fasta"
        self.DNA_DB_Split       = database_path + "ChocoPhlAn/ChocoPhlAn_split/"#"microbial_cds_db/microbial_all_cds.fasta"
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
        self.WEVOTEDB           = database_path + "WEVOTE_db/" #points to the location of taxdump, needed by WEVOTE
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
        self.Needle             = tool_path + "EMBOSS-6.6.0/emboss/needle"
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
        self.chart                      = script_path   + "pie_visualization.py"
        self.EC_Annotation_Post         = script_path   + "EC_Consolidation.py"

        self.Annotated_taxid            = script_path   + "Read_Classification/Get_TaxID.py"
        self.Constrain_classification   = script_path   + "Read_Classification/Constrain_Classification.py"
        self.Classification_combine     = script_path   + "Read_Classification/Combine_WEVOTE.py"
