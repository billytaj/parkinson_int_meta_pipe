FROM continuumio/anaconda3
MAINTAINER Billy Taj
RUN apt-get update
RUN apt-get install wget
RUN echo "hello from custom dockerfile"
WORKDIR /pipeline
ADD . /pipeline
CMD ["python", "hey.py"]

RUN wget https://github.com/weizhongli/cdhit/releases/download/V4.6.8/cd-hit-v4.6.8-2017-1208-source.tar.gz -O cdhit.tar.gz
RUN tar -xzvf cdhit.tar.gz

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.36.zip -O trimmomatic.zip

RUN apt-get install unzip
RUN unzip trimmomatic.zip 
RUN ls

