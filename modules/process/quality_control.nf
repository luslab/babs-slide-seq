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

// number of reads for each the barcodes matched or not, mapped or not
process reads_per_barcode_fastq {

	label "python"
	label "quality_control"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(fastq), path(script)

	output:
		tuple val(metadata), file("${name}.reads_per_barcode_fastq.csv")

	script:		

		name = metadata["name"]
		
		"""
		python3 $script $fastq "${name}.reads_per_barcode_fastq.csv"
		"""
}

