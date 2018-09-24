# Parkinson Lab Integrated Meta Transcriptomics/Genomics Pipeline
The Parkinson Lab Integrated Meta Transcriptomics/Genomics Pipeline is a software tool that will perform Transcriptomics and Genomic analysis
on Paired Reads of FASTQ data.

# How to install
This package is meant to work in conjunction with Docker/Singularity.  All of the prerequisite tools, including the pipeline code
is delivered via the Docker Hub. https://hub.docker.com/r/billyc59/parkinson_pipeline/
Alternatively, individual parts of the pipeline are avaiable from this Github repository.

# How to use
This pipeline comes with a config.ini file.  The user is meant to change, configure, and contort the file to point to the location of local files and Databases.
Our config file is written with Python's ConfigParser, using the basic interpretation.  
The following is an outline of wwhat each of the sections mean:

## Paramters
1) Output_Folder
Output_Folder is where the user would indicate to the pipeline where you want the output files to be dumped
2) Threads
Threads is the number of threads that the pipeline is allowed to use.  The pipeline is dependent on threads and parallelization to operate efficiently

## Sequences
1) Single
Single is for single reads.  Only fill this in if your sequence is a single file
2) Pair 1 and Pair 2
Pair 1 is for Forward Reads.  Pair 2 is for Reverse Reads.  Some portions of the code rely on the quality of the Forward Read file as a leading indicator of quality for the data filtering processes.  It is imperative that Pair 1 be used for your Forward Reads only, else you run the risk of a bad analysis. 

## Databases
We intentionally did not include any database files in this distribution so as to allow the user more flexibility in how they want to use the pipeline.  Below is a description of each of the fields we used.  
1) database_path
This field isn't apart of the parameters that the pipeline accepts.  It's a shortcut argument that makes filling the path to each database easier.
2) Univec_Core
The Univec_Core Database is used in the Vector Contaminents removal stage.  A copy can be found at: https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/
3) Adapter
The Adapter Database is used by Trimmomatic to trim Adapter segments out of the sequence files
A copy can be found inside the Trimmomatic tool installer, located at: http://www.usadellab.org/cms/?page=trimmomatic
This pipeline was built and tested using the TruSeq3-PE-2.fa Adapter Database
4) Host
The Host Database is used to filter out Host Contaminents from your seuqence file.  You will need to change this with the CDS database of whichever animal was used in your experiment.
We get our CDS databases from the NCBI, eg: ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human
5) Rfam
The Rfam Database is used by Infernal, the rRNA filter.  
A copy can be found here: http://rfam.xfam.org/
6) DNA_DB
The DNA DB is what we use to annotate the sequence data against.  We use the ChocoPhlAn database.
A copy can be found at: http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz
7) DNA_DB_Split
The ChocoPhlAn database is big.  We split it up and process the chunks simultaneously.  The pipeline will split it and dump the chunks at this location
8) Prot_DB
The Prot_DB is the protein db.  We use the non-redundant database from NCBI.  It will need to be indexed by DIAMOND before usage. (see DIAMOND for more details: https://github.com/bbuchfink/diamond)
It can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/db/
9) accession2taxid
This database links each accession to a taxid number.  It's used as part of a custom program in the pipeline.
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
10) nodes
This file is used in various parts in the pipeline.
It can be found at: ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
11) names
It can be found at the same location as nodes (10)
12) Kaiju_db
The Kaiju Database is used by Kaiju for taxonomic annotation.  Kaiju requires that the database be indexed before usage.  
Please see: https://github.com/bioinformatics-centre/kaiju/blob/master/README.md for more information.  Once the indexing is complete, the path to the database needs to be provided in this location
13) Centrifuge_db
The Centrifuge Database is used for the Centrifuge tool, which is a part of the enzyme annotation phase.  More information can be found here: https://ccb.jhu.edu/software/centrifuge/manual.shtml
We use the Nucleotide Database, after it has been indexed.  Details can be found at the link to Centrifuge

