import java.nio.file.Paths

process BAM_METRICS {
	label "python"
	label "bam_merics"
	label "process_low"
	tag { "${name}" }

	input:
	tuple val(metadata), path(bam), path(bai)

	output:
	tuple val(metadata), path("*.csv")

	script:		
	name = metadata["name"]
	def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
	"""
	reads_up_matching.py $bam "${name}.${suffix}.csv"
	"""
}

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

