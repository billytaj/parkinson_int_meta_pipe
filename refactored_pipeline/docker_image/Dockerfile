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
RUN mv bwa-0.7.17/* BWA
RUN rm -r bwa-0.7.17
RUN rm bwa.tar.bz2
WORKDIR /pipeline/BWA
RUN make
WORKDIR /pipeline

#------------------------------------------------------
# Install SAMTools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
RUN mkdir samtools
RUN mv samtools-1.7/* samtools
RUN rm -r samtools-1.7
RUN rm samtools.tar.bz2
WORKDIR /pipeline/samtools
RUN make
WORKDIR /pipeline

#----------------------------------------------------
# Install CD-HIT

ENV cdhit=cd-hit-v4.6.8-2017-1208
RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -O cdhit.tar.gz
RUN tar -xzvf cdhit.tar.gz
RUN mkdir cdhit
RUN mv cd-hit-v4.6.8-2017-1208/* cdhit/
RUN rm -r cd-hit-v4.6.8-2017-1208
RUN rm cdhit.tar.gz
WORKDIR /pipeline/cdhit
RUN make
WORKDIR /pipeline

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
