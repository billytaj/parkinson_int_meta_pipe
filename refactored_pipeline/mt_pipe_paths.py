
script_path = "/home/j/jparkins/billyc59/parkinson_int_meta_pipe/"

cdhit_dup = "/home/j/jparkins/mobolaji/Tools/CDHIT/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup"
Timmomatic = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
AdapterRemoval = "/home/j/jparkins/mobolaji/Tools/AdapterRemoval/adapterremoval-2.2.0/build/AdapterRemoval"
vsearch = "/home/j/jparkins/mobolaji/Tools/vsearch/vsearch-2.4.2-linux-x86_64/bin/vsearch"
UniVec_Core = "/home/j/jparkins/mobolaji/Databases/UniVec_Core.fasta"
Flash = "/home/j/jparkins/mobolaji/Tools/Flash/FLASH-1.2.11/flash"
Perl = "/home/j/jparkins/mobolaji/perl"
Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
Python = "/home/j/jparkins/mobolaji/python"
BWA = "/home/j/jparkins/mobolaji/Tools/BWA/bwa-0.7.5a/bwa"
SAMTOOLS = "/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/samtools"
BLAT = "/home/j/jparkins/mobolaji/Tools/pBLAT/pblat/pblat"
DIAMOND = "/home/j/jparkins/mobolaji/Tools/Diamond/diamond"
Blastdbcmd = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/blastdbcmd"
Makeblastdb = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/makeblastdb"
DNA_DB = "/scratch/j/jparkins/mobolaji/Microbial_cds_db/microbial_all_cds.fasta"
DNA_DB_Prefix = os.path.splitext(DNA_DB)[0]
DNA_DB_Extension = os.path.splitext(DNA_DB)[1]
Prot_DB = "/scratch/j/jparkins/mobolaji/NCBI_nr_db/nr"
Host = "/home/j/jparkins/mobolaji/Databases/Mouse_cds.fasta"
BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
Fastqc = "/home/j/jparkins/mobolaji/Tools/FastQC/fastqc"
Adapter = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmscan"
Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"
Filter_rRNA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
Reduplicate = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Reduplicate.py"
Map_reads_contigs = script_path + "Map_read_contigs.py"
Paired_Reads_Filter = script_path + "Paired_Reads_Filter.py"
BLAT_Contaminant_Filter = script_path + "BLAT_Contaminant_Filter.py"
File_splitter = script_path + "File_splitter.py"
Sort_Reads = script_path + "Read_Classification/Sort_Reads.py"
rRNA_Split_Jobs = script_path + "rRNA_Split_Jobs.py"
Map_reads_gene_BWA = script_path + "Map_read_gene_BWA.py"
Map_reads_gene_BLAT = script_path + "Map_read_gene_BLAT.py"
Map_reads_prot_DMND = script_path + "Map_read_prot_DMND.py"
Spades = "/home/j/jparkins/mobolaji/Tools/SPAdes/SPAdes-3.9.1-Linux/bin/spades.py"

EC_Annotation_Prep = script_path + "EC_Prediction_Scripts/0_Preprocess_Input.py"
Detect_Submit = script_path + "EC_Prediction_Scripts/1-1a_Detect_Submission.py"
EC_Annotation_Post = script_path + "EC_Prediction_Scripts/4a_EC_Consolidation.py"
Detect = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
Priam = "/home/j/jparkins/mobolaji/Tools/PRIAM/PRIAM_search.jar"
BLAST_dir = "/home/j/jparkins/mobolaji/Tools/BLAST/blast-2.2.26/bin/"
SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"

Nodes = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
Names = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
Annotated_taxid = script_path + "Read_Classification/Get_TaxID.py"
Contrain_classification = script_path + "Read_Classification/Constrain_Classification.py"
Classification_combine = script_path + "Read_Classification/Combine_WEVOTE.py"
WEVOTE = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/run_WEVOTE_PIPELINE.sh"
WEVOTEDB = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/WEVOTEDB"
accession2taxid = "/scratch/j/jparkins/mobolaji/accession2taxid/accession2taxid"
Kaiju = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju"
Kaiju2krona = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju2krona"
ktImportText = "/home/j/jparkins/mobolaji/Tools/Krona/KronaTools-2.7/bin/ktImportText"
Centrifuge = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge"
Centrifuge_report = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge-kreport"
kSLAM = "/home/j/jparkins/mobolaji/Tools/k-SLAM/k-SLAM/SLAM"

RPKM = script_path + "RPKM.py"