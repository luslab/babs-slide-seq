#!/bin/sh

bcl2fastq \
	--runfolder-dir $data_dir \
	--output-dir . \
	--sample-sheet $sample_sheet \
	--ignore-missing-bcls \
	--ignore-missing-filter \
	--barcode-mismatches $barcode_mismatches \
	--loading-threads $loading_threads \
	--processing-threads $processing_threads \
	--writing-threads $writing_threads

