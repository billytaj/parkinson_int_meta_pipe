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
# Install vsearch
RUN wget https://github.com/torognes/vsearch/archive/v2.7.1.tar.gz -O vsearch.tar.gz
RUN mkdir vsearch
RUN tar -xzvf vsearch.tar.gz
RUN mv vsearch-2.7.1/* vsearch
RUN rm -r vsearch-2.7.1
RUN rm vsearch.tar.gz
WORKDIR /pipeline/vsearch
RUN sh autogen.sh
RUN sh configure
RUN make
#RUN make -C vsearch
WORKDIR /pipeline

#------------------------------------------------------
# Install FLASH (Fast length adjustment of Short Reads)
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz -O FLASH.tar.gz
RUN tar -xzvf FLASH.tar.gz
RUN mkdir FLASH
RUN mv FLASH-1.2.11/* FLASH/
RUN rm -r FLASH-1.2.11
RUN rm FLASH.tar.gz
WORKDIR /pipeline/FLASH
RUN make
WORKDIR /pipeline

#------------------------------------------------------
# Install BWA
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 -O bwa.tar.bz2
RUN mkdir BWA
RUN tar -xjvf bwa.tar.bz2 
WORKDIR /pipeline/bwa-0.7.17
RUN make
WORKDIR /pipeline
RUN mv bwa-0.7.17/bwa BWA/
RUN rm -r bwa-0.7.17
RUN rm bwa.tar.bz2

#------------------------------------------------------
# Install SAMTools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
RUN mkdir samtools
#RUN mv samtools-1.7/* samtools
WORKDIR /pipeline/samtools-1.7
RUN make
RUN mv samtools /pipeline/samtools
WORKDIR /pipeline
RUN rm -r samtools-1.7
RUN rm samtools.tar.bz2

#----------------------------------------------------
# Install CD-HIT-DUP (from auxtools)

ENV cdhit=cd-hit-v4.6.8-2017-1208
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -O cdhit.tar.gz
RUN tar -xzvf cdhit.tar.gz
RUN mkdir cdhit_dup
WORKDIR /pipeline/cd-hit-v4.6.8-2017-1208/
RUN make
RUN mv cd-hit-auxtools/cd-hit-dup /pipeline/cdhit_dup/
WORKDIR /pipeline
RUN rm -r cd-hit-v4.6.8-2017-1208
RUN rm cdhit.tar.gz

#---------------------------------------------------
# Install trimmomatic

ENV trimm=trimmomatic-0.36
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip -O Trimmomatic.zip

RUN unzip Trimmomatic.zip 
RUN mkdir Trimmomatic
RUN mv Trimmomatic-0.36/* Trimmomatic
RUN rm -r Trimmomatic-0.36
RUN rm Trimmomatic.zip


#----------------------------------------------------
# Install AdapaterRemoval
RUN wget https://github.com/MikkelSchubert/adapterremoval/archive/v2.1.7.tar.gz -O adapterremoval.tar.gz
RUN mkdir adapterremoval
RUN tar -xzvf adapterremoval.tar.gz
RUN mv adapterremoval-2.1.7/* adapterremoval/
RUN rm -r adapterremoval-2.1.7
RUN rm adapterremoval.tar.gz
WORKDIR /pipeline/adapterremoval
RUN make 
WORKDIR /pipeline


#---------------------------------------------------
# Install BLAT (pblat)
RUN wget https://github.com/icebert/pblat/tarball/master -O PBLAT.tar.gz
RUN mkdir PBLAT
RUN tar -xzvf PBLAT.tar.gz
RUN mv icebert-pblat-e05e284/* PBLAT/
RUN rm -r icebert-pblat-e05e284
RUN rm PBLAT.tar.gz
WORKDIR /pipeline/PBLAT
RUN make
WORKDIR /pipeline

#-----------------------------------------------------
# Install DIAMOND
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.19/diamond-linux64.tar.gz -O DIAMOND.tar.gz
RUN mkdir DIAMOND
RUN tar -xzvf DIAMOND.tar.gz -C DIAMOND
RUN rm DIAMOND.tar.gz

#------------------------------------
# Install BLAST+
# compiling from source takes a much longer time than downloading the executables
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.7.1+-x64-linux.tar.gz -O BLAST_p.tar.gz
RUN tar -xzvf BLAST_p.tar.gz
RUN mkdir BLAST_p
RUN mv ncbi-blast-2.7.1+/bin/* BLAST_p
RUN rm BLAST_p.tar.gz
RUN rm -r ncbi-blast-2.7.1+

#-------------------------------------------
# Install Infernal
RUN wget http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz -O infernal.tar.gz
RUN tar -xzvf infernal.tar.gz
RUN mkdir infernal
RUN mv infernal-1.1.2-linux-intel-gcc/* infernal
RUN rm -r infernal-1.1.2-linux-intel-gcc
RUN rm infernal.tar.gz

#-------------------------------------------
# Install Kaiju
RUN wget https://github.com/bioinformatics-centre/kaiju/releases/download/v1.6.2/kaiju-1.6.2-linux-x86_64.tar.gz -O kaiju.tar.gz
RUN tar -xzvf kaiju.tar.gz
RUN mkdir kaiju
RUN mv kaiju-v1.6.2-linux-x86_64-static/bin/* kaiju/

#---------------------------------------
# Install KronaTools (for ktImportText)
# note, ktImportText is now a perl script, not a binary
RUN wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar -O KronaTools.tar
RUN tar -xvf KronaTools.tar
RUN echo "help"
RUN mkdir KronaTools
RUN mv KronaTools-2.7/* KronaTools/
WORKDIR KronaTools
RUN perl install.pl
WORKDIR /pipeline

#----------------------------------------
# Install k-SLAM
RUN wget https://github.com/aindj/k-SLAM/archive/master.zip -O k_slam.zip
RUN unzip k_slam.zip
RUN mkdir k-SLAM
RUN mv k-SLAM-master/* k-SLAM/
RUN rm k_slam.zip
RUN rm -r k-SLAM-master
WORKDIR k-SLAM 
RUN make -C build
RUN mv build/SLAM /pipeline
WORKDIR /pipeline
RUN rm -r k-SLAM/*
RUN mv SLAM k-SLAM/
WORKDIR /pipeline



#---------------------------------------
# Install Barrnap
RUN wget https://github.com/tseemann/barrnap/archive/master.zip -O barrnap.zip
RUN unzip barrnap.zip 
RUN mkdir barrnap
WORKDIR barrnap-master
RUN make 
RUN mv bin/barrnap /pipeline/barrnap
WORKDIR /pipeline
RUN rm -r barrnap-master
RUN rm barrnap.zip

#--------------------------------------------
# Install PRIAM
RUN mkdir PRIAM_search
RUN wget http://priam.prabi.fr/utilities/PRIAM_search.jar
RUN mv PRIAM_search.jar PRIAM_search

#------------------------------------
# Install WEVOTE
RUN wget https://github.com/aametwally/WEVOTE/archive/master.zip -O wevote.zip
RUN unzip wevote.zip
RUN mkdir WEVOTE
RUN mv WEVOTE-master/* WEVOTE/
RUN rm -r WEVOTE-master
RUN rm wevote.zip

