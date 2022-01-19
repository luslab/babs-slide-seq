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

//process tag_bam {
//
//	label "tagging"
//	label "sequencing"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam), path(script)
//
//	output:
//		tuple val(metadata), file("${name}.tagged.bam")
//
//	script:		
//		
//		name = metadata["name"]
//
//		"""
//		./$script $bam "${name}.tagged.bam"
//		"""
//}
//process primary {
//
//	label "tagging"
//	label "sequencing"
//
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam), path(script)
//
//	output:
//		tuple val(metadata), file("${name}.primary.bam")
//
//	script:		
//		
//		name = metadata["name"]
//
//		"""
//		./$script $bam "${name}.primary.bam"
//		"""
//}
//
//process multimap {
//
//	label "tagging"
//	label "sequencing"
//
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam), path(script)
//
//	output:
//		tuple val(metadata), file("${name}.multimap.bam")
//
//	script:		
//		
//		name = metadata["name"]
//
//		"""
//		./$script $bam "${name}.multimap.bam"
//		"""
//}
//
//process select {
//
//	label "tagging"
//	label "sequencing"
//
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam), path(script)
//
//	output:
//		tuple val(metadata), file("${name}.select.bam")
//
//	script:		
//		
//		name = metadata["name"]
//
//		"""
//		./$script $bam "${name}.select.bam"
//		"""
//}

