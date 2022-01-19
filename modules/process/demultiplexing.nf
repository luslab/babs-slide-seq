import java.nio.file.Paths

process bcl2fastq {

	//label "demultiplexing"

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if (filename.indexOf("Reports")!=-1 || filename.indexOf("Stats")!=-1)
			{
				"qc/demultiplexing/${filename}"
			}
			else
			{
				"files/demultiplexing/${filename}"
			}
		}

	input:
		tuple file(data_dir), path(sample_sheet)

	output:
		path "*.fastq.gz", emit: fastqs
		path "Reports", emit: reports
		path "Stats", emit: stats

	script:

		processing_threads = task.cpus
		loading_threads = 8
		writing_threads = 6
		barcode_mismatches = 0

		template "bcl2fastq.sh"
}

process merge_lanes {

	label "demultiplexing"

	tag { "${name}" }

	input:
		tuple val(metadata), path(fastq1), path(fastq2)

	output:
		tuple \
			val(metadata),
			path("${name}_R1.fastq.gz"),
			path("${name}_R2.fastq.gz")

	script:

		name = metadata["name"]

		"""
		zcat $fastq1 | gzip -c > "${name}_R1.fastq.gz"
		zcat $fastq2 | gzip -c > "${name}_R2.fastq.gz"
		"""
}

