import java.nio.file.Paths

process FASTQC {

	label "sequencing"
	//label "quality_control"
	tag { "${name}" }

	input:
		tuple val(metadata), path(fastq1), path(fastq2)

	output:
		tuple val(metadata), file("*_fastqc.html"), emit: html
		tuple val(metadata), file("*_fastqc.zip"), emit: zip

	script:		
		name = metadata["name"]
		
		"""
		fastqc $fastq1 $fastq2
		"""
}

process DUPLICATES {

	label "python"
	//label "quality_control"
	
	tag { "${basename}" }

	input:
		tuple val(metadata), path(csv), path(fastq1), path(fastq2), path(script)

	output:
		tuple val(metadata), file("${basename}.csv")

	script:		

		name = metadata["name"]
		status = metadata["status"]
		basename = "${name}.${status}"
		
		"""
		python3 $script $csv $fastq1 $fastq2 "${basename}.csv"
		"""
}

