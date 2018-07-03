FROM continuumio/anaconda3
MAINTAINER Billy Taj
RUN apt-get update
RUN apt-get install wget
RUN apt-get install unzip
RUN apt-get -y install g++
RUN apt-get -y install gcc
RUN apt-get -y install make
RUN apt-get -y install libz-dev
RUN apt-get -y install libbz2-dev
RUN apt-get -y install autoconf
RUN apt-get -y install build-essential
RUN apt-get -y install libboost-all-dev
RUN apt-get -y install flex
RUN apt-get -y install cmake
RUN apt-get -y install libncurses5-dev
RUN apt-get -y install liblzma-dev
RUN apt-get -y install perl
RUN conda update -n base conda -y
RUN conda install Biopython -y 
RUN conda install Pandas -y
RUN echo "hello from custom dockerfile"

WORKDIR /pipeline_tools


#----------------------------------------------------
# Install AdapaterRemoval
RUN wget https://github.com/MikkelSchubert/adapterremoval/archive/v2.1.7.tar.gz -O adapterremoval.tar.gz
RUN mkdir /pipeline_tools/adapterremoval
RUN tar -xzvf adapterremoval.tar.gz && rm adapterremoval.tar.gz
WORKDIR adapterremoval-2.1.7
RUN make 
RUN mv build/AdapterRemoval /pipeline_tools/adapterremoval/
WORKDIR /pipeline_tools
RUN rm -r adapterremoval-2.1.7 

#---------------------------------------
# Install Barrnap
RUN wget https://github.com/tseemann/barrnap/archive/master.zip -O barrnap.zip
RUN unzip barrnap.zip && rm barrnap.zip | mkdir barrnap
WORKDIR barrnap-master
RUN make | mv bin/barrnap /pipeline_tools/barrnap
WORKDIR /pipeline_tools
RUN rm -r barrnap-master

#------------------------------------
# Install BLAST+
# compiling from source takes a much longer time than downloading the executables
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz -O BLAST_p.tar.gz  
RUN tar --remove-files -xzvf BLAST_p.tar.gz && rm BLAST_p.tar.gz
RUN mkdir BLAST_p && mv ncbi-blast-2.7.1+/bin/* BLAST_p && rm -r ncbi-blast-2.7.1+

#---------------------------------------------------
# Install BLAT (pblat)
RUN wget https://github.com/icebert/pblat/tarball/master -O PBLAT.tar.gz | mkdir PBLAT
RUN tar --remove-files -xzvf PBLAT.tar.gz && rm PBLAT.tar.gz
WORKDIR /pipeline_tools/icebert-pblat-e05e284
RUN make && mv pblat /pipeline_tools/PBLAT/
WORKDIR /pipeline_tools
RUN rm -r icebert-pblat-e05e284

#------------------------------------------------------
# Install BWA
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 -O bwa.tar.bz2
RUN mkdir BWA && tar --remove-files -xjvf bwa.tar.bz2 && rm bwa.tar.bz2
WORKDIR /pipeline_tools/bwa-0.7.17
RUN make
WORKDIR /pipeline_tools
RUN mv bwa-0.7.17/bwa BWA/ && rm -r bwa-0.7.17


#----------------------------------------------------
# Install CD-HIT-DUP (from auxtools)
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -O cdhit.tar.gz
RUN tar --remove-files -xzvf cdhit.tar.gz && rm cdhit.tar.gz && mkdir cdhit_dup
WORKDIR cd-hit-v4.6.8-2017-1208/
RUN make && mv cd-hit-auxtools/cd-hit-dup /pipeline_tools/cdhit_dup/
WORKDIR /pipeline_tools
RUN rm -r cd-hit-v4.6.8-2017-1208



#----------------------------------------------------
# Install Centrifuge
RUN wget https://github.com/infphilo/centrifuge/archive/v1.0.3.tar.gz -O centrifuge.tar.gz
RUN tar --remove-files -xzvf centrifuge.tar.gz && rm centrifuge.tar.gz && mkdir /pipeline_tools/centrifuge
WORKDIR centrifuge-1.0.3
RUN make && mv centrifuge ../centrifuge && mv centrifuge-class ../centrifuge && mv centrifuge-kreport ../centrifuge && mv centrifuge-build ../centrifuge
WORKDIR /pipeline_tools
RUN rm -r centrifuge-1.0.3


#-----------------------------------------------------
# Install DIAMOND
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.19/diamond-linux64.tar.gz -O DIAMOND.tar.gz
RUN mkdir DIAMOND
RUN tar --remove-files -xzvf DIAMOND.tar.gz -C DIAMOND && rm DIAMOND.tar.gz


#------------------------------------------------------
# Install FLASH (Fast length adjustment of Short Reads)
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz -O FLASH.tar.gz
RUN tar --remove-files -xzvf FLASH.tar.gz && rm FLASH.tar.gz && mkdir FLASH
WORKDIR FLASH-1.2.11
RUN make && mv flash /pipeline_tools/FLASH
WORKDIR /pipeline_tools
RUN rm -r FLASH-1.2.11


#-------------------------------------------
# Install Infernal
RUN wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz -O infernal.tar.gz
RUN tar --remove-files -xzvf infernal.tar.gz && rm infernal.tar.gz && mkdir infernal
RUN mv infernal-1.1.2-linux-intel-gcc/binaries/* infernal && rm -r infernal-1.1.2-linux-intel-gcc
WORKDIR /pipeline_tools

#-------------------------------------------
# Install Kaiju
RUN wget https://github.com/bioinformatics-centre/kaiju/releases/download/v1.6.2/kaiju-1.6.2-linux-x86_64.tar.gz -O kaiju.tar.gz
RUN tar --remove-files -xzvf kaiju.tar.gz && rm kaiju.tar.gz
RUN mkdir kaiju && mv kaiju-v1.6.2-linux-x86_64-static/bin/* kaiju/ && rm -r kaiju-v1.6.2-linux-x86_64-static
WORKDIR /pipeline_tools


#---------------------------------------
# Install KronaTools (for ktImportText)
# note, ktImportText is now a perl script, not a binary.  We probably need most of this, if not all of this.  No cleaning
RUN wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar -O KronaTools.tar
RUN tar --remove-files -xvf KronaTools.tar && rm KronaTools.tar
RUN mkdir KronaTools && mv KronaTools-2.7/* KronaTools/ && rm -r KronaTools-2.7
WORKDIR KronaTools
RUN perl install.pl
WORKDIR /pipeline_tools

#----------------------------------------
# Install k-SLAM
RUN wget https://github.com/aindj/k-SLAM/archive/master.zip -O k_slam.zip
RUN unzip k_slam.zip && rm k_slam.zip && mkdir k-SLAM
WORKDIR k-SLAM-master
RUN make -C build && mv build/SLAM /pipeline_tools/k-SLAM/
WORKDIR /pipeline_tools
RUN rm -r k-SLAM-master

#--------------------------------------------
# Install PRIAM
RUN mkdir PRIAM_search
RUN wget http://priam.prabi.fr/utilities/PRIAM_search.jar
RUN mv PRIAM_search.jar PRIAM_search
WORKDIR /pipeline_tools


#------------------------------------------------------
# Install SAMTools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.8/samtools-1.8.tar.bz2 -O samtools.tar.bz2
RUN tar --remove-files -xjvf samtools.tar.bz2 && rm samtools.tar.bz2 && mkdir samtools
WORKDIR samtools-1.8
RUN make && mv samtools /pipeline_tools/samtools
WORKDIR /pipeline_tools
RUN rm -r samtools-1.8

#-------------------------------------------------------
# Install SPAdes
RUN wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz -O SPAdes.tar.gz
RUN tar --remove-files -xzvf SPAdes.tar.gz && rm SPAdes.tar.gz && mkdir SPAdes
RUN mv SPAdes-3.11.1-Linux/* SPAdes
RUN rm -r SPAdes-3.11.1-Linux


#---------------------------------------------------
# Install Trimmomatic and get the adapters files
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip -O Trimmomatic.zip
RUN unzip Trimmomatic.zip && rm Trimmomatic.zip && mkdir Trimmomatic
RUN mv Trimmomatic-0.36/* Trimmomatic && rm -r Trimmomatic-0.36


#----------------------------------------------------
# Install vsearch
RUN wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz -O vsearch.tar.gz
RUN mkdir vsearch && tar --remove-files -xzvf vsearch.tar.gz && rm vsearch.tar.gz 
WORKDIR vsearch-2.7.1
RUN sh autogen.sh && sh configure && make && mv bin/vsearch /pipeline_tools/vsearch
WORKDIR /pipeline_tools
RUN rm -r vsearch-2.7.1

#------------------------------------
# Install WEVOTE
RUN wget https://github.com/aametwally/WEVOTE/archive/master.zip -O wevote.zip
RUN unzip wevote.zip && rm wevote.zip
WORKDIR WEVOTE-master
RUN mkdir ../WEVOTE 
RUN ./install.sh && mkdir WEVOTE_db 
RUN mv bin/WEVOTE ../WEVOTE && mv downloadTaxonomyDB.sh ../WEVOTE && mkdir ../WEVOTE/WEVOTE_db && sh ../WEVOTE/downloadTaxonomyDB.sh ../WEVOTE/WEVOTE_db
WORKDIR /pipeline_tools
RUN rm -r WEVOTE-master

#--------------------------------------------
# get reference files
WORKDIR /pipeline_reference_files
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core -O UniVec_Core.fasta
RUN mkdir Trimmomatic_adapters
RUN mv /pipeline_tools/Trimmomatic/adapters/* /pipeline_reference_files/Trimmomatic_adapters
RUN wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz -O Mouse_cds.gz
RUN gunzip Mouse_cds.gz && mv Mouse_cds Mouse_cds.fasta
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/13.0/Rfam.cm.gz
RUN gunzip Rfam.cm.gz
RUN /pipeline_tools/infernal/cmpress Rfam.cm

WORKDIR /
CMD ["bash"]