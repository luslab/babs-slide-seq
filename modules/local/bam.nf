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

process BAM_FILTER {
	label "samtools"
	label 'process_low'
	
	tag { "${name}" }

	input:
	tuple val(metadata), path(bam), path(bai)

	output:
	tuple val(metadata), path("*.bam"), path("*.bam.bai")

	script:		
	name = metadata["name"]
	def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
	def expr    = task.ext.expr ?: ''
	"""
	echo "Filtering..."
	samtools view --expr '${expr}' --output "${name}.${suffix}.bam" $bam
	echo "Indexing..."
	samtools index "${name}.${suffix}.bam"
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

