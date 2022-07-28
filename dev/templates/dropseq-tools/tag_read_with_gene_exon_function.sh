#!/bin/sh

cp -v $refflat annot.refFlat

drop-seq \
	TagReadWithGeneExonFunction \
	INPUT=$bam \
	OUTPUT=$out_bam \
	TAG=GE \
	ALLOW_MULTI_GENE_READS=false \
	ANNOTATIONS_FILE=annot.refFlat

