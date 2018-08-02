# Parkinson Lab Integrated Meta Transcriptomics/Genomics Pipeline
The Parkinson Lab Integrated Meta Transcriptomics/Genomics Pipeline is a software tool that will perform Transcriptomics and Genomic analysis
on Paired Reads of FASTQ data.

# How to install
This package is meant to work in conjunction with a Singularity Machine.  All of the prerequisite tools, including the pipeline code
is delivered to you via the Docker Hub. https://hub.docker.com/r/billyc59/parkinson_pipeline/
Alternatively, you can pull individual parts from this Github repository.

Note:  there are issues with the pipeline if you are not a privileged user.  The python scripts will try to write.  If it is not allowed, you
will need to copy the entire contents of /raw_pipeline to a location where you are free to write.

