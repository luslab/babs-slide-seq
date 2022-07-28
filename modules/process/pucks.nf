process SHUFFLING {
	tag { "${name}" }

	label "pucks"
	label "python"
	label "ultra_low"

	input:
	tuple val(name), val(barcodes)

	output:
	tuple val(name), val("ordered"), path("${name}.ordered.txt"), emit: ordered
	tuple val(name), val("shuffled"), path("${name}.shuffled.txt"), emit: shuffled

	script:
	"""
	shuffling.py "${barcodes}" "${name}"
	"""
}
