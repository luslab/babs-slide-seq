import java.nio.file.Paths

process MARK_DUPLICATES {
	//label "tagging"
	label "sequencing"
	label 'process_high'

	tag { "${name}" }

	publishDir Paths.get( params.outdir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf(".txt") != -1 )
			{
				"qc/${filename}"
			}
			else
			{
				"files/${filename}"
			}
		}

	input:
		tuple val(metadata), path(bam)
	
	output:
		tuple val(metadata), path("${name}.dup.bam"), emit: bam
		tuple val(metadata), path("${name}.dup.txt"), emit: metrics
	
	script:
		
		name = metadata["name"]

		"""
		picard-tools \
			MarkDuplicates \
				I=$bam \
				O=${name}.dup.bam \
				M=${name}.dup.txt
		"""
}

process BAM_TAG {

	label "tagging"
	label "sequencing"
	
	tag { "${basename}" }

	input:
		tuple val(metadata), path(bam), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"


		"""
		./$script $bam "${basename}.bam"
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process BAM_TAH_HMEM {

	label "tagging"
	label "sequencing"
	
	tag { "${basename}" }

	input:
		tuple val(metadata), path(bam), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"


		"""
		./$script $bam "${basename}.bam"
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process umis_per_barcode {

	label "tagging"
	label "sequencing"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), path(bai), path(script)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai")

	script:		
		
		name = metadata["name"]
		suffix = "umis"
		basename = "${name}.${suffix}"
		threshold = params.umis_threshold

		"""
		./$script --threshold $threshold $bam "${basename}.bam"
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

process htseq {

	label "tagging"
	label "sequencing"

	tag { "${name}" }
	
	input:
		tuple val(metadata), path(bam), path(bai)

	output:
		tuple val(metadata), path("${out_bam}"), path("${out_bam}.bai"), emit: bam
		tuple val(metadata), path("${out_txt}"), emit: txt

	script:		
		
		name = metadata["name"]
		out_sam = "${name}.htseq.sam"
		out_bam = "${name}.htseq.bam"
		out_txt = "${name}.htseq.txt"
		gtf = metadata["gtf"]

		"""
		htseq-count --format bam --samout sample.sam $bam $gtf > $out_txt
		samtools view -H $bam > $out_sam
		cat sample.sam >> $out_sam
		samtools view -S -b $out_sam > $out_bam
		samtools index $out_bam
		"""
}

process select {

	label "tagging"
	label "sequencing"
	
	tag { "${basename}" }

	publishDir Paths.get( params.out_dir ),
		mode: "copy",
		overwrite: "true",
		saveAs: { filename ->
			if ( filename.indexOf(".csv") != -1 )
			{
				"qc/${filename}"
			}
			else
			{
				"files/${filename}"
			}
		}

	input:
		tuple val(metadata), path(bam), val(suffix), path(script)

	output:
		tuple val(metadata), path("${basename}.bam"), path("${basename}.bam.bai"), emit: bam
		tuple val(metadata), val("unique"), path("${basename}.unique.csv"), emit: unique_reads
		tuple val(metadata), val("resolved"), path("${basename}.resolved.csv"), emit: resolved_reads
		tuple val(metadata), val("unresolved"), path("${basename}.unresolved.csv"), emit: unresolved_reads

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"


		"""
		./$script "${basename}" $bam "${basename}.bam"
		echo "Indexing..."
		samtools index "${basename}.bam"
		"""
}

