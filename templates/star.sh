#!/bin/sh

STAR \
	--runMode alignReads \
	--runThreadN ${task.cpus} \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesCommand zcat \
	--genomeDir $genome_index \
	--outFileNamePrefix ${name}. \
	--readFilesIn $fastq

