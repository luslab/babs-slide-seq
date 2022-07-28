process PLOT {
	label 'plot'
	tag 'FILL_TAG_IN'

	input:
	tuple val(metadata), val(value)

	output:
	tuple val(metadata), file("*.png"), emit: png
	tuple val(metadata), file("*.pdf"), emit: pdf

	script:
	name = metadata["name"]

	def args    = task.ext.args ?: ''
	def script    = task.ext.script ?: ''
    def ext     = task.ext.ext ?: 'txt'
    """
    python3 $script $args $value ${name}.${ext}
    """
}


// process plot_1_arg {

// 	label "plot"
	
// 	tag { "${name}" }

// 	input:
// 		tuple val(metadata), path(arg), val(suffix), path(script)

// 	output:
// 		tuple val(metadata), file("${basename}.png"), emit: png
// 		tuple val(metadata), file("${basename}.pdf"), emit: pdf

// 	script:		
		
// 		name = metadata["name"]
// 		basename = "${name}.${suffix}"

// 		"""
// 		python3 $script $arg $basename
// 		"""
// }

// process plot_1_arg_1_val {

// 	label "plot"
	
// 	tag { "${name}" }

// 	input:
// 		tuple val(metadata), path(arg), val(value), val(suffix), path(script)

// 	output:
// 		tuple val(metadata), file("${basename}.png"), emit: png
// 		tuple val(metadata), file("${basename}.pdf"), emit: pdf

// 	script:		
		
// 		name = metadata["name"]
// 		basename = "${name}.${suffix}"

// 		"""
// 		python3 $script $arg $value $basename
// 		"""
// }

// process plot_2_args {

// 	label "plot"
	
// 	tag { "${name}" }

// 	input:
// 		tuple val(metadata), path(arg1), path(arg2), val(suffix), path(script)

// 	output:
// 		tuple val(metadata), file("${basename}.png"), emit: png
// 		tuple val(metadata), file("${basename}.pdf"), emit: pdf

// 	script:		
		
// 		name = metadata["name"]
// 		basename = "${name}.${suffix}"

// 		"""
// 		python3 $script $arg1 $arg2 $basename
// 		"""
// }

// process plot_1_val {

// 	label "plot"
	
// 	tag { "${name}" }

// 	input:
// 		tuple val(metadata), path(arg), val(value), val(suffix), path(script)

// 	output:
// 		tuple val(metadata), file("${basename}.png"), emit: png
// 		tuple val(metadata), file("${basename}.pdf"), emit: pdf

// 	script:		
		
// 		name = metadata["name"]
// 		basename = "${name}.${suffix}"

// 		"""
// 		python3 $script $arg $basename $value
// 		"""
// }

