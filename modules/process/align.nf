import java.nio.file.Paths

process star {

	tag { "${name}" }
	time "10:00:00"

	publishDir Paths.get( params.out_dir ),
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
	
	output:
		tuple val(metadata), path("${name}.Log.*"), emit: log
		tuple val(metadata), path("${name}.SJ.out.tab"), emit: tab
		tuple \
			val(metadata),
			path("${name}.Aligned.sortedByCoord.out.bam"), emit: bam

	script:		
		
		name = metadata["name"]
		index = metadata["genome"]

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

