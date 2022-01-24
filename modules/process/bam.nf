import java.nio.file.Paths

//process bam_index {
//
//	label "sequencing"
//	label "bam"
//	
//	tag { "${name}" }
//
//	publishDir Paths.get( params.out_dir , "files" ),
//		mode: "copy",
//		overwrite: "true"
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), path("${new_bam}"), path("${new_bam}.bai")
//
//	script:		
//		
//		name = metadata["name"]
//		new_bam = bam.getFileName().toString().replace(".bam", "") + ".index.bam"
//
//		"""
//		cp -v $bam $new_bam
//		samtools index $new_bam
//		"""
//}

process bam_filter {

	label "samtools"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "files" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), val(expr)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		samtools view --expr '${expr}' --output "${basename}.bam" $bam
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process bam_metrics {

	label "python"
	label "bam_merics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "qc" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		"""
}

process bam_metrics_hmem {

	label "python"
	label "bam_merics"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir , "qc" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(bam), path(bai), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.csv")

	script:		

		name = metadata["name"]
		basename = "${name}.${suffix}"
		
		"""
		python3 $script $bam "${basename}.csv"
		"""
}

