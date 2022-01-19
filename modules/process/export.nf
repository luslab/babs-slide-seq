
process merge_plots {

	label "export"
	label "python_2"
	
	tag { "${name}" }

	input:
		tuple val(name), path(pdfs)

	output:
		file "${name}.pdf"

	script:		
		"""
		pdfunite $pdfs "${name}.pdf"
		"""
}

process rename_coords {

	label "export"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(csv)

	output:
		file "${name}.csv"

	script:		
		
		name = metadata["name"]

		"""
		cp -v $csv "${name}.csv"
		"""
}

process dge {

	label "export"
	label "sequencing"

	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), path(script)

	output:
		tuple val(metadata), file("${directory}")

	script:		
		
		name = metadata["name"]
		gtf = metadata["gtf"]
		directory = "${name}_dge"

		"""
		./$script $gtf $bam .

		mkdir $directory
		cat matrix.mtx | gzip -c > $directory/matrix.mtx.gz
		cat features.tsv | gzip -c > $directory/features.tsv.gz
		cat barcodes.tsv | gzip -c > $directory/barcodes.tsv.gz
		"""
}

