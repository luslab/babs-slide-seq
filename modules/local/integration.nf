process GET_BARCODES {
	label "integration"
	label "process_low"
	tag { "${name}" }
	container 'ubuntu:20.04'

	input:
	tuple val(metadata), path(csv)

	output:
	tuple val(metadata), file("${name}.barcodes.txt")

	script:			
	name = metadata["name"]
	"""
	cat $csv \
		| sed 's/,.*//g' \
		| sort \
		| uniq \
		> "${name}.barcodes.txt"
	"""
}

process HAMMING {
	label "integration"
	label "process_low"
	tag { "${name}" }
	container 'chrischeshire/slideseq-hamming:latest'

	input:
	tuple val(metadata), path(read_barcodes), path(puck_barcodes)

	output:
	tuple val(metadata), file("*.csv")

	script:			
	name = metadata["name"]
	bcd = metadata["barcodes"]
	csv = "${name}.${bcd}.hamming.csv"
	"""
	hostname
	nvidia-smi
	hamming $read_barcodes $puck_barcodes "${csv}"
	"""
}

process MATCHER {
	label "integration"
	label "python"
	label "process_low"
	tag { "${name}" }

	input:
	tuple val(metadata), path(hamming), path(reads), path(puck)

	output:
	tuple val(metadata), file("*.map.matching.csv"    ), emit: mapping
	tuple val(metadata), file("*values.matching.csv"  ), emit: values
	tuple val(metadata), file("*.metrics.matching.csv"), emit: metrics
	tuple val(metadata), file("*.csv"                 ), emit: coords

	script:				
	name = metadata["name"]
	shuf = metadata["barcodes"]
	basename = name + "." + shuf
	"""
	matcher.py $hamming $reads $puck "${basename}" \
		$params.barcode_errors_threshold \
		$params.barcode_max_matches \
		$params.barcode_max_entropy
	"""
}

process ADD_MATCH {
	label "integration"
	label "sequencing"
	label "process_low"
	tag { "${name}" }

	input:
		tuple val(metadata), path(mapping), path(bam)

	output:
	tuple val(metadata), file("${name}.matched.bam"), file("${name}.matched.bam.bai")

	script:		
	name = metadata["name"]

	"""
	./$script $mapping $bam "${name}.matched.bam"
	samtools index "${name}.matched.bam"
	"""
}
