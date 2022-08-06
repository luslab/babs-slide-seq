
process MERGE_PLOTS {
	label "export"
	label "python"
	
	tag { "${name}" }

	input:
	tuple val(name), path(pdfs)

	output:
	file "*.pdf"

	script:		
	"""
	pdfunite $pdfs "${name}.pdf"
	"""
}

process RENAME_COORDS {
	label "export"	
	tag { "${name}" }

	input:
	tuple val(metadata), path(csv)

	output:
	file "*.csv"

	script:		
	name = metadata["name"]
	"""
	cp -v $csv "${name}.csv"
	"""
}

process DGE {
	label "export"
	label "sequencing"

	tag { "${name}" }

	input:
	tuple val(metadata), path(bam)
	path gtf

	output:
	tuple val(metadata), file("${directory}")

	script:
	name = metadata["name"]
	directory = "${name}_dge"
	"""
	count --directory . $gtf $bam

	mkdir $directory
	cat matrix.mtx | gzip -c > $directory/matrix.mtx.gz
	cat features.tsv | gzip -c > $directory/features.tsv.gz
	cat barcodes.tsv | gzip -c > $directory/barcodes.tsv.gz
	"""
}

