import os.path
from configparser import ConfigParser, ExtendedInterpolation

class tool_path_obj:
    def __init__ (self, config_path):

        if config_path:
            config = ConfigParser() #change this to ex
            config.read(config_path)
        else:
            print("no config found, defaulting")
            config = None

        script_path             = "/pipeline/Scripts"
        tool_path               = "/pipeline_tools/"
        database_path           = "/project/j/jparkin/Lab_Databases/"
        custom_database_path    = "/pipeline/custom_databases/"

        #----------------------------------------------------------
        # Reference Databases
        # Note: default host is Mouse CDS
        
        if config:
            self.UniVec_Core        = config["Databases"]["UniVec_Core"]        if config["Databases"]["UniVec_Core"]       or config["Databases"]["UniVec_Core"]       == "" else os.path.join(database_path, "univec_core/UniVec_Core.fasta")
            self.Adapter            = config["Databases"]["Adapter"]            if config["Databases"]["Adapter"]           or config["Databases"]["Adapter"]           == "" else os.path.join(database_path, "Trimmomatic_adapters/TruSeq3-PE-2.fa")
            self.Host               = config["Databases"]["Host"]               if config["Databases"]["Host"]              or config["Databases"]["Host"]              == "" else os.path.join(database_path, "Mouse_cds/Mouse_cds.fasta")
            self.Rfam               = config["Databases"]["Rfam"]               if config["Databases"]["Rfam"]              or config["Databases"]["Rfam"]              == "" else os.path.join(database_path, "Rfam/Rfam.cm")
            self.DNA_DB             = config["Databases"]["DNA_DB"]             if config["Databases"]["DNA_DB"]            or config["Databases"]["DNA_DB"]            == "" else os.path.join(database_path, "ChocoPhlAn/ChocoPhlAn.fasta")
            self.DNA_DB_Split       = config["Databases"]["DNA_DB_Split"]       if config["Databases"]["DNA_DB_Split"]      or config["Databases"]["DNA_DB_Split"]      == "" else os.path.join(database_path, "ChocoPhlAn/ChocoPhlAn_split/")
            self.Prot_DB            = config["Databases"]["Prot_DB"]            if config["Databases"]["Prot_DB"]           or config["Databases"]["Prot_DB"]           == "" else os.path.join(database_path, "nr/nr")
            self.Prot_DB_reads      = config["Databases"]["Prot_DB_reads"]            if config["Databases"]["Prot_DB_reads"]           or config["Databases"]["Prot_DB_reads"]           == "" else os.path.join(database_path, "nr/nr")
            self.accession2taxid    = config["Databases"]["accession2taxid"]    if config["Databases"]["accession2taxid"]   or config["Databases"]["accession2taxid"]   == "" else os.path.join(database_path, "accession2taxid/accession2taxid")
            self.nodes              = config["Databases"]["nodes"]              if config["Databases"]["nodes"]             or config["Databases"]["nodes"]             == "" else os.path.join(database_path, "WEVOTE_db/nodes.dmp")
            self.names              = config["Databases"]["names"]              if config["Databases"]["names"]             or config["Databases"]["names"]             == "" else os.path.join(database_path, "WEVOTE_db/names.dmp")
            self.Kaiju_db           = config["Databases"]["Kaiju_db"]           if config["Databases"]["Kaiju_db"]          or config["Databases"]["Kaiju_db"]          == "" else os.path.join(database_path, "kaiju_db/kaiju_db_nr.fmi")
            self.Centrifuge_db      = config["Databases"]["Centrifuge_db"]      if config["Databases"]["Centrifuge_db"]     or config["Databases"]["Centrifuge_db"]     == "" else os.path.join(database_path, "centrifuge_db/nt")
            self.SWISS_PROT         = config["Databases"]["SWISS_PROT"]         if config["Databases"]["SWISS_PROT"]        or config["Databases"]["SWISS_PROT"]        == "" else os.path.join(database_path, "swiss_prot_db/swiss_prot_db")
            self.SWISS_PROT_map     = config["Databases"]["SWISS_PROT_map"]     if config["Databases"]["SWISS_PROT_map"]    or config["Databases"]["SWISS_PROT_map"]    == "" else os.path.join(database_path, "swiss_prot_db/SwissProt_EC_Mapping.tsv")
            self.PriamDB            = config["Databases"]["PriamDB"]            if config["Databases"]["PriamDB"]           or config["Databases"]["PriamDB"]           == "" else os.path.join(database_path, "PRIAM_db/")
            self.DetectDB           = config["Databases"]["DetectDB"]           if config["Databases"]["DetectDB"]          or config["Databases"]["DetectDB"]          == "" else os.path.join(database_path, "DETECTv2")
            self.WEVOTEDB           = config["Databases"]["WEVOTEDB"]           if config["Databases"]["WEVOTEDB"]          or config["Databases"]["WEVOTEDB"]          == "" else os.path.join(database_path, "WEVOTE_db/")
            self.EC_pathway         = config["Databases"]["EC_pathway"]         if config["Databases"]["EC_pathway"]        or config["Databases"]["EC_pathway"]        == "" else os.path.join(database_path, "EC_pathway.txt")
            self.path_to_superpath  = config["Databases"]["path_to_superpath"]  if config["Databases"]["path_to_superpath"] or config["Databases"]["path_to_superpath"] == "" else os.path.join(custom_database_path, "pathway_to_superpathway.csv")
            self.mgm_model          = config["Databases"]["MetaGeneMark_model"] if config["Databases"]["MetaGeneMark_model"]or config["Databases"]["MetaGeneMark_model"]== "" else os.path.join(tool_path, "mgm/MetaGeneMark_v1.mod")
            
            
        else:
            self.UniVec_Core        = os.path.join(database_path,           "univec_core/UniVec_Core.fasta")
            self.Adapter            = os.path.join(database_path,           "Trimmomatic_adapters/TruSeq3-PE-2.fa")
            self.Host               = os.path.join(database_path,           "Human_cds/Human_cds.fasta")
            self.Rfam               = os.path.join(database_path,           "Rfam/Rfam.cm")
            self.DNA_DB             = os.path.join(database_path,           "ChocoPhlAn/ChocoPhlAn.fasta")
            self.DNA_DB_Split       = os.path.join(database_path,           "ChocoPhlAn/ChocoPhlAn_split/")
            self.Prot_DB            = os.path.join(database_path,           "nr/nr")
            self.Prot_DB_reads      = os.path.join(database_path,           "nr/nr")
            self.accession2taxid    = os.path.join(database_path,           "accession2taxid/accession2taxid")
            self.nodes              = os.path.join(database_path,           "WEVOTE_db/nodes.dmp")
            self.names              = os.path.join(database_path,           "WEVOTE_db/names.dmp")
            self.Kaiju_db           = os.path.join(database_path,           "kaiju_db/kaiju_db_nr.fmi")
            self.Centrifuge_db      = os.path.join(database_path,           "centrifuge_db/nt")
            self.SWISS_PROT         = os.path.join(database_path,           "swiss_prot_db/swiss_prot_db")
            self.SWISS_PROT_map     = os.path.join(database_path,           "swiss_prot_db/SwissProt_EC_Mapping.tsv")
            self.PriamDB            = os.path.join(database_path,           "PRIAM_db/")
            self.DetectDB           = os.path.join(database_path,           "DETECTv2")
            self.WEVOTEDB           = os.path.join(database_path,           "WEVOTE_db/")
            self.EC_pathway         = os.path.join(custom_database_path,    "EC_pathway.txt")
            self.path_to_superpath  = os.path.join(custom_database_path,    "pathway_to_superpathway.csv")
            self.mgm_model          = os.path.join(tool_path,               "mgm/MetaGeneMark_v1.mod")
   
        try:
            self.enzyme_db = config["Databases"]["enzyme_db"] if config["Databases"]["enzyme_db"] or config["Databases"]["enzyme_db"] == "" else os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt")
        except:
            self.enzyme_db = os.path.join(custom_database_path, "FREQ_EC_pairs_3_mai_2020.txt")
            
        #----------------------------------------------------------
        # external tools
        
        if config:
            self.Python         = config["Tools"]["Python"]         if config["Tools"]["Python"]            or config["Tools"]["Python"]            == "" else "python3"
            self.Java           = config["Tools"]["Java"]           if config["Tools"]["Java"]              or config["Tools"]["Java"]              == "" else "java -Xmx10000m -jar"
            self.cdhit_dup      = config["Tools"]["cdhit_dup"]      if config["Tools"]["cdhit_dup"]         or config["Tools"]["cdhit_dup"]         == "" else os.path.join(tool_path, "cdhit_dup/cd-hit-dup")
            self.Timmomatic     = config["Tools"]["Timmomatic"]     if config["Tools"]["Timmomatic"]        or config["Tools"]["Timmomatic"]        == "" else os.path.join(tool_path, "Trimmomatic/trimmomatic-0.36.jar")
            self.AdapterRemoval = config["Tools"]["AdapterRemoval"] if config["Tools"]["AdapterRemoval"]    or config["Tools"]["AdapterRemoval"]    == "" else os.path.join(tool_path, "adapterremoval/AdapterRemoval")
            self.vsearch        = config["Tools"]["vsearch"]        if config["Tools"]["vsearch"]           or config["Tools"]["vsearch"]           == "" else os.path.join(tool_path, "vsearch/vsearch")
            self.Flash          = config["Tools"]["Flash"]          if config["Tools"]["Flash"]             or config["Tools"]["Flash"]             == "" else os.path.join(tool_path, "FLASH/flash")
            self.BWA            = config["Tools"]["BWA"]            if config["Tools"]["BWA"]               or config["Tools"]["BWA"]               == "" else os.path.join(tool_path, "BWA/bwa")
            self.SAMTOOLS       = config["Tools"]["SAMTOOLS"]       if config["Tools"]["SAMTOOLS"]          or config["Tools"]["SAMTOOLS"]          == "" else os.path.join(tool_path, "samtools/samtools")
            self.BLAT           = config["Tools"]["BLAT"]           if config["Tools"]["BLAT"]              or config["Tools"]["BLAT"]              == "" else os.path.join(tool_path, "PBLAT/pblat")
            self.DIAMOND        = config["Tools"]["DIAMOND"]        if config["Tools"]["DIAMOND"]           or config["Tools"]["DIAMOND"]           == "" else os.path.join(tool_path, "DIAMOND/diamond")
            self.Blastp         = config["Tools"]["Blastp"]         if config["Tools"]["Blastp"]            or config["Tools"]["Blastp"]            == "" else os.path.join(tool_path, "BLAST_p/blastp")
            self.Needle         = config["Tools"]["Needle"]         if config["Tools"]["Needle"]            or config["Tools"]["Needle"]            == "" else os.path.join(tool_path, "EMBOSS-6.6.0/emboss/stretcher")
            self.Blastdbcmd     = config["Tools"]["Blastdbcmd"]     if config["Tools"]["Blastdbcmd"]        or config["Tools"]["Blastdbcmd"]        == "" else os.path.join(tool_path, "BLAST_p/blastdbcmd")
            self.Makeblastdb    = config["Tools"]["Makeblastdb"]    if config["Tools"]["Makeblastdb"]       or config["Tools"]["Makeblastdb"]       == "" else os.path.join(tool_path, "BLAST_p/makeblastdb")
            self.Barrnap        = config["Tools"]["Barrnap"]        if config["Tools"]["Barrnap"]           or config["Tools"]["Barrnap"]           == "" else os.path.join(tool_path, "Barrnap/bin/barrnap")
            self.Infernal       = config["Tools"]["Infernal"]       if config["Tools"]["Infernal"]          or config["Tools"]["Infernal"]          == "" else os.path.join(tool_path, "infernal/cmsearch")
            self.Kaiju          = config["Tools"]["Kaiju"]          if config["Tools"]["Kaiju"]             or config["Tools"]["Kaiju"]             == "" else os.path.join(tool_path, "kaiju/kaiju")
            self.Centrifuge     = config["Tools"]["Centrifuge"]     if config["Tools"]["Centrifuge"]        or config["Tools"]["Centrifuge"]        == "" else os.path.join(tool_path, "centrifuge/centrifuge")
            self.Priam          = config["Tools"]["Priam"]          if config["Tools"]["Priam"]             or config["Tools"]["Priam"]             == "" else os.path.join(tool_path, "PRIAM_search/PRIAM_search.jar")
            self.Detect         = config["Tools"]["Detect"]         if config["Tools"]["Detect"]            or config["Tools"]["Detect"]            == "" else os.path.join(script_path, "Detect_2.2.7.py")
            self.BLAST_dir      = config["Tools"]["BLAST_dir"]      if config["Tools"]["BLAST_dir"]         or config["Tools"]["BLAST_dir"]         == "" else os.path.join(tool_path, "BLAST_p")
            self.WEVOTE         = config["Tools"]["WEVOTE"]         if config["Tools"]["WEVOTE"]            or config["Tools"]["WEVOTE"]            == "" else os.path.join(tool_path, "WEVOTE/WEVOTE")
            self.Spades         = config["Tools"]["Spades"]         if config["Tools"]["Spades"]            or config["Tools"]["Spades"]            == "" else os.path.join(tool_path, "SPAdes/bin/spades.py")
            self.MetaGeneMark   = config["Tools"]["MetaGeneMark"]   if config["Tools"]["MetaGeneMark"]      or config["Tools"]["MetaGeneMark"]      == "" else os.path.join(tool_path, "mgm/gmhmmp")
        else:
            self.Python         = "python3"
            self.Java           = "java -jar"
            self.cdhit_dup      = os.path.join(tool_path, "cdhit_dup/cd-hit-dup")
            self.Timmomatic     = os.path.join(tool_path, "Trimmomatic/trimmomatic-0.36.jar")
            self.AdapterRemoval = os.path.join(tool_path, "adapterremoval/AdapterRemoval")
            self.vsearch        = os.path.join(tool_path, "vsearch/vsearch")
            self.Flash          = os.path.join(tool_path, "FLASH/flash")
            self.BWA            = os.path.join(tool_path, "BWA/bwa")
            self.SAMTOOLS       = os.path.join(tool_path, "samtools/samtools")
            self.BLAT           = os.path.join(tool_path, "PBLAT/pblat")
            self.DIAMOND        = os.path.join(tool_path, "DIAMOND/diamond")
            self.Blastp         = os.path.join(tool_path, "BLAST_p/blastp")
            self.Needle         = os.path.join(tool_path, "EMBOSS-6.6.0/emboss/stretcher")
            self.Blastdbcmd     = os.path.join(tool_path, "BLAST_p/blastdbcmd")
            self.Makeblastdb    = os.path.join(tool_path, "BLAST_p/makeblastdb")
            self.Barrnap        = os.path.join(tool_path, "Barrnap/bin/barrnap")
            self.Infernal       = os.path.join(tool_path, "infernal/cmsearch")
            self.Kaiju          = os.path.join(tool_path, "kaiju/kaiju")
            self.Centrifuge     = os.path.join(tool_path, "centrifuge/centrifuge")
            self.Priam          = os.path.join(tool_path, "PRIAM_search/PRIAM_search.jar")
            self.Detect         = os.path.join(script_path, "Detect_2.2.8.py")
            self.BLAST_dir      = os.path.join(tool_path, "BLAST_p")
            self.WEVOTE         = os.path.join(tool_path, "WEVOTE/WEVOTE")
            self.Spades         = os.path.join(tool_path, "SPAdes/bin/spades.py")
            self.MetaGeneMark   = os.path.join(tool_path, "mgm/gmhmmp")
            
        #--------------------------------------------
        # Python scripts

        self.sam_trimmer                = os.path.join(script_path, "read_sam.py")
        self.sort_reads                 = os.path.join(script_path, "read_sort.py")
        self.duplicate_repopulate       = os.path.join(script_path, "read_repopulation.py")
        self.orphaned_read_filter       = os.path.join(script_path, "read_orphan.py")
        self.remove_tag                 = os.path.join(script_path, "read_remove_tag.py")
        self.BLAT_Contaminant_Filter    = os.path.join(script_path, "read_BLAT_filter.py")
        self.File_splitter              = os.path.join(script_path, "read_split.py")
        self.barrnap_post               = os.path.join(script_path, "read_rRNA_barrnap.py")
        self.rRNA_filter                = os.path.join(script_path, "read_rRNA_infernal.py")
        self.map_contig                 = os.path.join(script_path, "assembly_map.py")
        self.contig_duplicate_remover   = os.path.join(script_path, "assembly_deduplicate.py")
        self.Map_reads_gene_BWA         = os.path.join(script_path, "ga_BWA_generic.py")
        self.Map_reads_gene_BLAT        = os.path.join(script_path, "ga_BLAT_generic.py")
        self.Map_reads_prot_DMND        = os.path.join(script_path, "ga_Diamond_generic.py")
        self.GA_final_merge             = os.path.join(script_path, "ga_Final_merge.py")
        self.EC_Annotation_Post         = os.path.join(script_path, "ea_combine_v4.py")
        self.Annotated_taxid            = os.path.join(script_path, "ta_taxid.py")
        self.Constrain_classification   = os.path.join(script_path, "ta_constrain_taxonomy_v2.py")
        self.Classification_combine     = os.path.join(script_path, "ta_combine_v2.py")
        self.Wevote_parser              = os.path.join(script_path, "ta_wevote_parser.py")
        self.taxa_table                 = os.path.join(script_path, "output_taxa_table.py")
        self.RPKM                       = os.path.join(script_path, "output_table_v3.py")
        self.format_RPKM                = os.path.join(script_path, "output_reformat_rpkm_table.py")
        self.read_count                 = os.path.join(script_path, "output_read_counts.py")
        self.read_quality_metrics       = os.path.join(script_path, "output_read_quality_metrics.py")
        self.contig_stats               = os.path.join(script_path, "output_contig_stats.py")
        self.ec_heatmap                 = os.path.join(script_path, "output_EC_metrics.py")
        self.data_change_metrics        = os.path.join(script_path, "output_data_change_metrics.py")
        self.get_unique_host_reads      = os.path.join(script_path, "output_get_host_reads.py")
        self.remove_gaps_in_fasta       = os.path.join(script_path, "remove_gaps_in_fasta.py")
        self.parse_sam                  = os.path.join(script_path, "output_parse_sam.py")
        self.are_you_in_a_contig        = os.path.join(script_path, "output_are_you_in_a_contig.py")
        self.convert_contig_segments    = os.path.join(script_path, "output_convert_gene_map_contig_segments.py")
        self.output_filter_taxa         = os.path.join(script_path, "output_filter_taxa.py")
        self.output_filter_ECs          = os.path.join(script_path, "output_filter_ECs.py")
        #--------------------------------------------------
        # miscellaneous values
        BWA_mem_default = 50
        BLAT_mem_default = 10 #100MB
        DIAMOND_mem_default = 50 #60GB
        DETECT_mem_default = 50
        Infernal_mem_default = 50
        
        chunk_size_default = 10000
        BWA_job_limit_default = 80
        BLAT_job_limit_default = 80
        DIAMOND_job_limit_default = 80
        DETECT_job_limit_default = 1000
        Infernal_job_limit_default = 1000
        
        keep_all_default = "yes"
        keep_quality_default = "no"
        keep_host_default = "no"
        keep_repop_default = "no"
        keep_vector_default = "no"
        keep_assemble_contigs_default = "no"
        keep_rRNA_default = "no"
        keep_GA_BWA_default = "no"
        keep_GA_BLAT_default = "no"
        keep_GA_DIAMOND_default = "no"
        keep_GA_final_default = "no"
        keep_TA_default = "no"
        keep_EC_default = "no"
        keep_outputs_default = "no"
        
        delete_all_default = "no"
        
        rRNA_job_delay_default = 5
        BWA_job_delay_default = 5
        BLAT_job_delay_default = 5
        DIAMOND_job_delay_default = 5
        DETECT_job_delay_default = 5
        
        if config:
            self.target_rank                = config["Settings"]["Target_Rank"]                 if config["Settings"]["Target_Rank"]                or config["Settings"]["Target_Rank"]                == "" else "genus"
            self.adapterremoval_minlength   = config["Settings"]["AdapterRemoval_minlength"]    if config["Settings"]["AdapterRemoval_minlength"]   or config["Settings"]["AdapterRemoval_minlength"]   == "" else 30
            self.show_unclassified          = config["Settings"]["Show_unclassified"]           if config["Settings"]["Show_unclassified"]          or config["Settings"]["Show_unclassified"]          == "" else "No"
            self.rpkm_cutoff                = config["Settings"]["RPKM_cutoff"]                 if config["Settings"]["RPKM_cutoff"]                or config["Settings"]["RPKM_cutoff"]                == "" else 0.01
            try:
                self.BWA_mem_threshold = config["Settings"]["BWA_mem_threshold"] if config["Settings"]["BWA_mem_threshold"] or config["Settings"]["BWA_mem_threshold"]  == "" else BWA_mem_default
            except KeyError:
                self.BWA_mem_threshold = BWA_mem_default
                
            try:
                self.BLAT_mem_threshold = config["Settings"]["BLAT_mem_threshold"] if config["Settings"]["BLAT_mem_threshold"] or config["Settings"]["BLAT_mem_threshold"] == "" else BLAT_mem_default
            except KeyError:
                self.BLAT_mem_threshold = BLAT_mem_default  
                
            try:
                self.DIAMOND_mem_threshold  = config["Settings"]["DIAMOND_mem_threshold"] if config["Settings"]["DIAMOND_mem_threshold"] or config["Settings"]["DIAMOND_mem_threshold"] == "" else DIAMOND_mem_default
            except KeyError:
                self.DIAMOND_mem_threshold = DIAMOND_mem_default
                
            try:
                self.DETECT_mem_threshold = config["Settings"]["DETECT_mem_threshold"] if config["Settings"]["DETECT_mem_threshold"] or config["Settings"]["DETECT_mem_threshold"] == "" else DETECT_mem_default
            except KeyError:
                self.DETECT_mem_threshold = DETECT_mem_default
            
            try:
                self.Infernal_mem_threshold = config["Settings"]["Infernal_mem_threshold"] if config["Settings"]["Infernal_mem_threshold"] or config["Settings"]["Infernal_mem_threshold"] == "" else Infernal_mem_default
            except KeyError:
                self.Infernal_mem_threshold = Infernal_mem_default
                
                
            try:
                self.chunk_size  = config["Settings"]["rRNA_chunk_size"] if config["Settings"]["rRNA_chunk_size"] or config["Settings"]["rRNA_chunk_size"] == "" else chunk_size_default
            except KeyError:
                self.chunk_size = chunk_size_default
                
            try:
                self.BWA_job_limit  = config["Settings"]["BWA_job_limit"] if config["Settings"]["BWA_job_limit"] or config["Settings"]["BWA_job_limit"] == "" else BWA_job_limit_default
            except KeyError:
                self.BWA_job_limit = BWA_job_limit_default
                
            try:
                self.BLAT_job_limit  = config["Settings"]["BLAT_job_limit"] if config["Settings"]["BLAT_job_limit"] or config["Settings"]["BLAT_job_limit"] == "" else BLAT_job_limit_default
            except KeyError:
                self.BLAT_job_limit = BLAT_job_limit_default
            
            try:
                self.DIAMOND_job_limit  = config["Settings"]["DIAMOND_job_limit"] if config["Settings"]["DIAMOND_job_limit"] or config["Settings"]["DIAMOND_job_limit"] == "" else DIAMOND_job_limit_default
            except KeyError:
                self.DIAMOND_job_limit = DIAMOND_job_limit_default
                
            try:
                self.DETECT_job_limit = config["Settings"]["DETECT_job_limit"] if config["Settings"]["DETECT_job_limit"] or config["Settings"]["DETECT_job_limit"] == "" else DETECT_job_limit_default
            except KeyError:
                self.DETECT_job_limit = DETECT_job_limit_default
            
            try:
                self.Infernal_job_limit = config["Settings"]["rRNA_job_limit"] if config["Settings"]["rRNA_job_limit"] or config["Settings"]["rRNA_job_limit"] == "" else Infernal_job_limit_default
            except KeyError:
                self.Infernal_job_limit = Infernal_job_limit_default
            
            #------------------------------------------------------------------------------------------------
            try:
                self.keep_all = config["Settings"]["keep_all"] if config["Settings"]["keep_all"] or config["Settings"]["keep_all"] == "" else keep_all_default
            except KeyError:
                self.keep_all = keep_all_default
                
            try:
                self.keep_quality = config["Settings"]["keep_quality"] if config["Settings"]["keep_quality"] or config["Settings"]["keep_quality"] == "" else keep_quality_default
            except KeyError:
                self.keep_quality = keep_quality_default
                
            try:
                self.keep_host = config["Settings"]["keep_host"] if config["Settings"]["keep_host"] or config["Settings"]["keep_host"] == "" else keep_host_default
            except KeyError:
                self.keep_host = keep_host_default
                    
            try:
                self.keep_vector = config["Settings"]["keep_vector"] if config["Settings"]["keep_vector"] or config["Settings"]["keep_vector"] == "" else keep_vector_default
            except KeyError:
                self.keep_vector = keep_vector_default
            
            try:
                self.keep_rRNA = config["Settings"]["keep_rRNA"] if config["Settings"]["keep_rRNA"] or config["Settings"]["keep_rRNA"] == "" else keep_rRNA_default
            except KeyError:
                self.keep_rRNA = keep_rRNA_default
            
            try:
                self.keep_repop = config["Settings"]["keep_repop"] if config["Settings"]["keep_repop"] or config["Settings"]["keep_repop"] == "" else keep_repop_default
            except KeyError:
                self.keep_repop = keep_repop_default
            
            try:
                self.keep_assemble_contigs = config["Settings"]["keep_assemble_contigs"] if config["Settings"]["keep_assemble_contigs"] or config["Settings"]["keep_assemble_contigs"] == "" else keep_assemble_contigs_default
            except KeyError:
                self.keep_assemble_contigs = keep_assemble_contigs_default
            
            try:
                self.keep_GA_BWA = config["Settings"]["keep_GA_BWA"] if config["Settings"]["keep_GA_BWA"] or config["Settings"]["keep_GA_BWA"] == "" else keep_GA_BWA_default
            except KeyError:
                self.keep_GA_BWA = keep_GA_BWA_default
            
            try:
                self.keep_GA_BLAT = config["Settings"]["keep_GA_BLAT"] if config["Settings"]["keep_GA_BLAT"] or config["Settings"]["keep_GA_BLAT"] == "" else keep_GA_BLAT_default
            except KeyError:
                self.keep_GA_BLAT = keep_GA_BLAT_default
            
            try:
                self.keep_GA_DIAMOND = config["Settings"]["keep_GA_DIAMOND"] if config["Settings"]["keep_GA_DIAMOND"] or config["Settings"]["keep_GA_DIAMOND"] == "" else keep_GA_DIAMOND_default
            except KeyError:
                self.keep_GA_DIAMOND = keep_GA_DIAMOND_default
                
            try:
                self.keep_GA_final = config["Settings"]["keep_GA_final"] if config["Settings"]["keep_GA_final"] or config["Settings"]["keep_GA_final"] == "" else keep_gene_annotation_FINAL_MERGE_default
            except KeyError:
                self.keep_GA_final = keep_GA_final_default    
                
            try:
                self.keep_TA = config["Settings"]["keep_taxonomic_annotation"] if config["Settings"]["keep_taxonomic_annotation"] or config["Settings"]["keep_taxonomic_annotation"] == "" else keep_TA_default
            except KeyError:
                self.keep_TA = keep_TA_default    
            
            try:
                self.keep_EC = config["Settings"]["keep_enzyme_annotation"] if config["Settings"]["keep_enzyme_annotation"] or config["Settings"]["keep_enzyme_annotation"] == "" else keep_EC_default
            except KeyError:
                self.keep_EC = keep_EC_default
             
            try:
                self.keep_outputs = config["Settings"]["keep_outputs"] if config["Settings"]["keep_outputs"] or config["Settings"]["keep_outputs"] == "" else keep_outputs_default
            except KeyError:
                self.keep_outputs = keep_outputs_default
            
            #------------------------------------------------------------------------
            try:
                self.rRNA_job_delay = config["Settings"]["rRNA_job_delay"] if config["Settings"]["rRNA_job_delay"] or config["Settings"]["rRNA_job_delay"] == "" else rRNA_job_delay_default
            except KeyError:
                self.rRNA_job_delay = rRNA_job_delay_default
                
            try:
                self.BWA_job_delay = config["Settings"]["BWA_job_delay"] if config["Settings"]["BWA_job_delay"] or config["Settings"]["BWA_job_delay"] == "" else BWA_job_delay_default
            except KeyError:
                self.BWA_job_delay = BWA_job_delay_default
                
            try:
                self.BLAT_job_delay = config["Settings"]["BLAT_job_delay"] if config["Settings"]["BLAT_job_delay"] or config["Settings"]["BLAT_job_delay"] == "" else BLAT_job_delay_default
            except KeyError:
                self.BLAT_job_delay = BLAT_job_delay_default
                
            try:
                self.DIAMOND_job_delay = config["Settings"]["DIAMOND_job_delay"] if config["Settings"]["DIAMOND_job_delay"] or config["Settings"]["DIAMOND_job_delay"] == "" else DIAMOND_job_delay_default
            except KeyError:
                self.DIAMOND_job_delay = DIAMOND_job_delay_default
                
            try:
                self.DETECT_job_delay = config["Settings"]["DETECT_job_delay"] if config["Settings"]["DETECT_job_delay"] or config["Settings"]["DETECT_job_delay"] == "" else DETECT_job_delay
            except KeyError:
                self.DETECT_job_delay = DETECT_job_delay_default
                
                
                
                
             
        else:
            self.target_rank = "genus"
            self.adapterremoval_minlength = 30
            self.show_unclassified = "No"
            self.rpkm_cutoff = 0.01
            self.BWA_mem_threshold = BWA_mem_default
            self.BLAT_mem_threshold = BLAT_mem_default
            self.DIAMOND_mem_threshold = DIAMOND_mem_default
            self.DETECT_mem_threshold = DETECT_mem_default
            self.Infernal_mem_threshold = Infernal_mem_default
            self.chunk_size = chunk_size_default
            self.BWA_job_limit = BWA_job_limit_default
            self.BLAT_job_limit = BLAT_job_limit_default
            self.DIAMOND_job_limit = DIAMOND_job_limit_default
            self.DETECT_job_limit = DETECT_job_limit_default
            self.Infernal_job_limit = Infernal_job_limit_default
            
            self.keep_all = keep_all_default
            self.keep_quality = keep_quality_default
            self.keep_host = keep_host_default
            self.keep_vector = keep_vector_default
            self.keep_rRNA = keep_rRNA_default
            self.keep_repop = self.keep_repop_default
            self.keep_assemble_contigs = keep_assemble_contigs_default
            self.keep_GA_BWA = keep_GA_BWA_default
            self.keep_GA_BLAT = keep_GA_BLAT_default
            self.keep_GA_DIAMOND = keep_GA_DIAMOND_default
            self.keep_GA_final = GA_final_default
            self.keep_TA = keep_TA_default
            self.keep_EC = keep_EC_default
            self.keep_outputs = keep_outputs_default
            
            self.rRNA_job_delay = rRNA_job_delay_default
            self.BWA_job_delay = BWA_job_delay_default
            self.BLAT_job_delay = BLAT_job_delay_default
            self.DIAMOND_job_delay = DIAMOND_job_delay_default
            self.DETECT_job_delay = DETECT_job_delay_default
            