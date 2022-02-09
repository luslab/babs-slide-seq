import java.nio.file.Paths

process fastqc {

	label "sequencing"
	//label "quality_control"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "qc" , "fastqc" ),
		mode: "copy",
		overwrite: "true"

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

