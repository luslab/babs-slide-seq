import java.nio.file.Paths

//process get_reads_per_barcode {
//
//	label "samtools"
//	label "integration"
//	
//	tag { "${name}" }
//
//	publishDir Paths.get( params.out_dir , "qc" ),
//		mode: "copy",
//		overwrite: "true"
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${csv}")
//
//	script:		
//		
//		name = metadata["name"]
//		csv = "${name}.reads_per_barcode.csv"
//
//		"""
//		samtools view \
//			--expr '[us] == "MATCHED" && [as] == "MAPPED" && [bc]' $bam \
//			| awk '{
//				v=\$0; gsub(".*bc:Z:", "", v);
//				gsub("\\t.*", "", v);
//				printf "%s\\t%s\\n", v, \$1 }
//			' \
//			| sort \
//			| uniq \
//			| awk '{ print \$1 }' \
//			| sort \
//			| uniq -c \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			| sort -t , -k 2,2rn \
//			> "${csv}"
//		"""
//}

//process get_reads_per_barcode {
//
//	label "python_2"
//	label "integration"
//	
//	tag { "${name}" }
//
//	publishDir Paths.get( params.out_dir , "qc" ),
//		mode: "copy",
//		overwrite: "true"
//
//	input:
//		tuple val(metadata), path(bam), path(bai), path(script)
//
//	output:
//		tuple val(metadata), file("${csv}")
//
//	script:		
//		
//		name = metadata["name"]
//		csv = "${name}.reads_per_barcode.csv"
//
//		"""
//		python3 $script $bam $csv
//		"""
//}

process get_barcodes {

	label "samtools"
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
		tuple val(metadata), file("${name}.matched.bam")

	script:		
		
		name = metadata["name"]

		"""
		./$script $mapping $bam "${name}.matched.bam"
		"""
}
