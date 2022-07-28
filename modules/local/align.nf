import java.nio.file.Paths

process STAR {
	tag { "${name}" }
	label 'process_high'
	container = 'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:afaaa4c6f5b308b4b6aa2dd8e99e1466b2a6b0cd-0'

	publishDir Paths.get( params.outdir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf(".Log.") != -1 )
			{
				"qc/align/${filename}"
			}
			else
			{
				"files/align/${filename}"
			}
		}

	input:
	tuple val(metadata), path(fastq)
	path(index)
	
	output:
	tuple val(metadata), path("${name}.Log.*"), emit: log
	tuple val(metadata), path("${name}.SJ.out.tab"), emit: tab
	tuple val(metadata), path("${name}.Aligned.sortedByCoord.out.bam"), emit: bam

	script:		
	name = metadata["name"]

	"""
	STAR \
		--runMode alignReads \
		--runThreadN $task.cpus \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outSAMattributes NH HI AS nM NM MD \
		--genomeDir $index \
		--outFileNamePrefix "${name}." \
		--readFilesCommand zcat \
		--readFilesIn $fastq
	"""
}

