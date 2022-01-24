import java.nio.file.Paths

process get_barcodes {

	label "integration"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "qc" ),
		mode: "copy",
		overwrite: "true"

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

process hamming {

	label "integration"
	label "gpu"
	
	time "06:00:00"
	memory "40G"

	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "qc" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple \
			val(metadata),
			path(read_barcodes),
			path(puck_barcodes),
			path(script), path(cl)

	output:
		tuple val(metadata), file("${csv}")

	script:		
		
		name = metadata["name"]
		bcd = metadata["barcodes"]
		basename = "${name}.${bcd}"
		csv = "${basename}.hamming.csv"

		"""
		./$script $read_barcodes $puck_barcodes "${csv}"
		"""
}

process matcher {

	label "integration"
	label "python"
	time "02:00:00"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "qc" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(hamming), path(reads), path(puck), path(script)

	output:
		tuple val(metadata), file("${basename}.map.matching.csv"), emit: mapping
		tuple val(metadata), file("${basename}.values.matching.csv"), emit: values
		tuple val(metadata), file("${basename}.metrics.matching.csv"), emit: metrics
		tuple val(metadata), file("${basename}.csv"), emit: coords

	script:		
		
		name = metadata["name"]
		shuf = metadata["barcodes"]
		basename = name + "." + shuf

		"""
		python3 $script $hamming $reads $puck "${basename}" \
			$params.barcode_errors_threshold \
			$params.barcode_max_matches \
			$params.barcode_max_entropy
		"""
}

process add_match {

	label "integration"
	label "sequencing"
	time "03:00:00"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "files" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(mapping), path(bam), path(script)

	output:
		tuple \
			val(metadata),
			file("${name}.matched.bam"),
			file("${name}.matched.bam.bai")

	script:		
		
		name = metadata["name"]

		"""
		./$script $mapping $bam "${name}.matched.bam"
		samtools index "${name}.matched.bam"
		"""
}
