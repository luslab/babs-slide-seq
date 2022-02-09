FROM ubuntu:21.04

RUN apt-get update
RUN apt-get --yes upgrade
RUN apt-get install --yes apt-utils

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London
RUN apt-get install --yes tzdata

RUN apt-get install --yes less

RUN apt-get install --yes \
	bowtie \
	bowtie-examples \
	bowtie2 \
	bowtie2-examples \
	bwa \
	clustalo \
	clustalw \
	clustalx \
	cutadapt \
	fastqc \
	igv \
	kallisto \
	kallisto-examples \
	libpicard-java \
	libpicard-java-doc \
	libseqan2-dev \
	libseqan3-dev \
	libarmadillo-dev \
	libarmadillo10 \
	libboost-all-dev \
	multiqc \
	picard \
	picard-tools \
	python-htseq \
	python3-cutadapt \
	python3-htseq \
	rna-star \
	rsem \
	salmon \
	samtools \
	seqan-apps \
	seqan3-doc \
	seqtk \
	trimmomatic

