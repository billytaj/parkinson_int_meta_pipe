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
RUN echo "hello from custom dockerfile"

WORKDIR /pipeline
ADD . /pipeline
CMD ["python", "hey.py"]

#----------------------------------------------------
# Install AdapaterRemoval
RUN wget https://github.com/MikkelSchubert/adapterremoval/archive/v2.1.7.tar.gz -O adapterremoval.tar.gz
RUN mkdir adapterremoval
RUN tar --remove-files -xzvf adapterremoval.tar.gz
WORKDIR /pipeline/adapterremoval-2.1.7
RUN make 
RUN mv build/AdapterRemoval /pipeline/adapterremoval/
WORKDIR /pipeline


#---------------------------------------
# Install Barrnap
RUN wget https://github.com/tseemann/barrnap/archive/master.zip -O barrnap.zip
RUN unzip barrnap.zip && rm barrnap.zip
RUN mkdir barrnap
WORKDIR barrnap-master
RUN make 
RUN mv bin/barrnap /pipeline/barrnap
WORKDIR /pipeline
RUN rm -r barrnap-master

#------------------------------------
# Install BLAST+
# compiling from source takes a much longer time than downloading the executables
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz -O BLAST_p.tar.gz
RUN tar --remove-files -xzvf BLAST_p.tar.gz
RUN mkdir BLAST_p
RUN mv ncbi-blast-2.7.1+/bin/* BLAST_p
RUN rm -r ncbi-blast-2.7.1+

#---------------------------------------------------
# Install BLAT (pblat)
RUN wget https://github.com/icebert/pblat/tarball/master -O PBLAT.tar.gz
RUN mkdir PBLAT
RUN tar --remove-files -xzvf PBLAT.tar.gz
WORKDIR /pipeline/icebert-pblat-e05e284
RUN make
RUN mv pblat /pipeline/PBLAT/
WORKDIR /pipeline
RUN rm -r icebert-pblat-e05e284

#------------------------------------------------------
# Install BWA
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 -O bwa.tar.bz2
RUN mkdir BWA
RUN tar --remove-files -xjvf bwa.tar.bz2 
WORKDIR /pipeline/bwa-0.7.17
RUN make
WORKDIR /pipeline
RUN mv bwa-0.7.17/bwa BWA/
RUN rm -r bwa-0.7.17

#----------------------------------------------------
# Install CD-HIT-DUP (from auxtools)
ENV cdhit=cd-hit-v4.6.8-2017-1208
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -O cdhit.tar.gz
RUN tar --remove-files -xzvf cdhit.tar.gz
RUN mkdir cdhit_dup
WORKDIR /pipeline/cd-hit-v4.6.8-2017-1208/
RUN make
RUN mv cd-hit-auxtools/cd-hit-dup /pipeline/cdhit_dup/
WORKDIR /pipeline
RUN rm -r cd-hit-v4.6.8-2017-1208

#-----------------------------------------------------
# Install DIAMOND
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.19/diamond-linux64.tar.gz -O DIAMOND.tar.gz
RUN mkdir DIAMOND
RUN tar --remove-files -xzvf DIAMOND.tar.gz -C DIAMOND

#------------------------------------------------------
# Install FLASH (Fast length adjustment of Short Reads)
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz -O FLASH.tar.gz
RUN tar --remove-files -xzvf FLASH.tar.gz
RUN mkdir FLASH
WORKDIR FLASH-1.2.11
RUN make
RUN mv flash /pipeline/FLASH
WORKDIR /pipeline
RUN rm -r FLASH-1.2.11

#-------------------------------------------
# Install Infernal
RUN wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz -O infernal.tar.gz
RUN tar --remove-files -xzvf infernal.tar.gz
RUN mkdir infernal
RUN mv infernal-1.1.2-linux-intel-gcc/* infernal
RUN rm -r infernal-1.1.2-linux-intel-gcc

#-------------------------------------------
# Install Kaiju
RUN wget https://github.com/bioinformatics-centre/kaiju/releases/download/v1.6.2/kaiju-1.6.2-linux-x86_64.tar.gz -O kaiju.tar.gz
RUN tar --remove-files -xzvf kaiju.tar.gz
RUN mkdir kaiju
RUN mv kaiju-v1.6.2-linux-x86_64-static/bin/* kaiju/

#---------------------------------------
# Install KronaTools (for ktImportText)
# note, ktImportText is now a perl script, not a binary
RUN wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar -O KronaTools.tar
RUN tar --remove-files -xvf KronaTools.tar
RUN echo "help"
RUN mkdir KronaTools
RUN mv KronaTools-2.7/* KronaTools/
WORKDIR KronaTools
RUN perl install.pl
WORKDIR /pipeline

#----------------------------------------
# Install k-SLAM
RUN wget https://github.com/aindj/k-SLAM/archive/master.zip -O k_slam.zip
RUN unzip k_slam.zip && rm k_slam.zip
RUN mkdir k-SLAM
WORKDIR k-SLAM-master
RUN make -C build
RUN mv build/SLAM /pipeline/k-SLAM/
WORKDIR /pipeline
RUN rm -r k-SLAM-master
#--------------------------------------------
# Install PRIAM
RUN mkdir PRIAM_search
RUN wget http://priam.prabi.fr/utilities/PRIAM_search.jar
RUN mv PRIAM_search.jar PRIAM_search

#------------------------------------------------------
# Install SAMTools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2 -O samtools.tar.bz2
RUN tar --remove-files -xjvf samtools.tar.bz2
RUN mkdir samtools
#RUN mv samtools-1.7/* samtools
WORKDIR /pipeline/samtools-1.7
RUN make
RUN mv samtools /pipeline/samtools
WORKDIR /pipeline
RUN rm -r samtools-1.7

#---------------------------------------------------
# Install Trimmomatic
ENV trimm=trimmomatic-0.36
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip -O Trimmomatic.zip
RUN unzip Trimmomatic.zip && rm Trimmomatic.zip
RUN mkdir Trimmomatic
RUN mv Trimmomatic-0.36/* Trimmomatic
RUN rm -r Trimmomatic-0.36

#----------------------------------------------------
# Install vsearch
RUN wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz -O vsearch.tar.gz
RUN mkdir vsearch
RUN tar --remove-files -xzvf vsearch.tar.gz
WORKDIR vsearch-2.7.1
RUN sh autogen.sh
RUN sh configure
RUN make
RUN mv bin/vsearch /pipeline/vsearch
WORKDIR /pipeline
RUN rm -r vsearch-2.7.1

#------------------------------------
# Install WEVOTE
RUN wget https://github.com/aametwally/WEVOTE/archive/master.zip -O wevote.zip
RUN unzip wevote.zip && rm wevote.zip
RUN mkdir WEVOTE
RUN mv WEVOTE-master/* WEVOTE/
RUN rm -r WEVOTE-master

