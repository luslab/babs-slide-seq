import java.nio.file.Paths

process mark_duplicates {

	//label "tagging"
	label "sequencing"
	time "10:00:00"

	tag { "${name}" }

	publishDir Paths.get( params.out_dir ),
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

process bam_tag {

	label "tagging"
	label "sequencing"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), val(suffix), path(script)

	output:
		tuple val(metadata), file("${name}.${suffix}.bam")

	script:		
		
		name = metadata["name"]


		"""
		./$script $bam "${name}.${suffix}.bam"
		"""
}

process bam_tag_hmem {

	label "tagging"
	label "sequencing"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), val(suffix), path(script)

	output:
		tuple val(metadata), file("${name}.${suffix}.bam")

	script:		
		
		name = metadata["name"]


		"""
		./$script $bam "${name}.${suffix}.bam"
		"""
}

process umis_per_barcode {

	label "tagging"
	label "sequencing"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), path(script)

	output:
		tuple val(metadata), file("${name}.umis.bam")

	script:		
		
		name = metadata["name"]
		threshold = params.umis_threshold


		"""
		./$script $bam "${name}.umis.bam" $threshold
		"""
}

process gene {

	label "tagging"
	label "sequencing"

	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), path(script)

	output:
		tuple val(metadata), file("${name}.gene.bam")

	script:		
		
		name = metadata["name"]
		gtf = metadata["gtf"]

		"""
		./$script $gtf $bam "${name}.gene.bam"
		"""
}

process dropseq {

	label "tagging"
	label "sequencing"

	tag { "${name}" }
	
	input:
		tuple val(metadata), path(bam)

	output:
		tuple val(metadata), path("${out_bam}")

	script:		
		
		name = metadata["name"]
		out_bam = "${name}.dropseq.bam"
		refflat = metadata["refflat"]

		template "dropseq-tools/tag_read_with_gene_exon_function.sh"
}

process dropseq_tag {

	label "tagging"
	label "sequencing"

	tag { "${name}" }

	input:
		tuple val(metadata), path(bam), path(script)

	output:
		tuple val(metadata), file("${name}.dropseqtag.bam")

	script:		
		
		name = metadata["name"]
		gtf = metadata["gtf"]

		"""
		./$script $gtf $bam "${name}.dropseqtag.bam"
		"""
}

