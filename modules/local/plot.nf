process PLOT {
	label 'plot'
	label 'ultra_low'
	tag "$name"

	input:
	tuple val(metadata), path(file), val(value)

	output:
	tuple val(metadata), file("*.png"), emit: png
	tuple val(metadata), file("*.pdf"), emit: pdf

	script:
	name = metadata["name"]

	def args    = task.ext.args ?: ''
	def script  = task.ext.script ?: ''
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    """
    $script $args $file $value ${name}.${suffix}
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

