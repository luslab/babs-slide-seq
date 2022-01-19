
process shuffling {

	tag { "${name}" }

	label "pucks"
	label "python"

	input:
		tuple val(name), path(barcodes), path(script)

	output:
		tuple \
			val(name),
			val("ordered"), path("${name}.ordered.txt"), emit: ordered
		tuple \
			val(name),
			val("shuffled"),
			path("${name}.shuffled.txt"), emit: shuffled

	script:		
		"""
		python3 $script $barcodes $name
		"""
}

