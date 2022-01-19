
process plot_1_arg {

	label "plot"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(arg), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg $basename
		"""
}

process plot_2_args {

	label "plot"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(arg1), path(arg2), val(suffix), path(script)

	output:
		tuple val(metadata), file("${basename}.png"), emit: png
		tuple val(metadata), file("${basename}.pdf"), emit: pdf

	script:		
		
		name = metadata["name"]
		basename = "${name}.${suffix}"

		"""
		python3 $script $arg1 $arg2 $basename
		"""
}

//process plot_histo_hamming {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(ordered_csv), path(shuffled_csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".histo_hamming"
//
//		"""
//		python3 $script $ordered_csv $shuffled_csv $basename
//		"""
//}
//
//process plot_spatial_umis {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(mtx), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".spatial_umis"
//
//		"""
//		python3 $script $mtx $csv $basename
//		"""
//}

//process plot_balance_barcode {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".balance_barcode"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_balance_umi {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".balance_umi"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_barcode_align {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".barcode_align"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_barcode_extraction {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".barcode_extraction"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_barcode_matching {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".barcode_matching"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_duplicates {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".duplicates"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_histo_errors {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".histo_errors"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_histo_function {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".histo_function"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_histo_genes {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(mtx), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".histo_genes"
//
//		"""
//		python3 $script $mtx $basename
//		"""
//}
//
//process plot_histo_umis {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(mtx), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".histo_umis"
//
//		"""
//		python3 $script $mtx $basename
//		"""
//}
//
//process plot_mappings {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".mappings"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_reads_fraction {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".reads_fraction"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_resolved {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".resolved"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_select {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".selected"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_umis_per_barcode {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(mtx), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".umis_per_barcode"
//
//		"""
//		python3 $script $mtx $basename
//		"""
//}
//
//process plot_up_align {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".up_align"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}
//
//process plot_up_matching {
//
//	label "plot"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(csv), path(script)
//
//	output:
//		tuple val(metadata), file("${basename}.png"), emit: png
//		tuple val(metadata), file("${basename}.pdf"), emit: pdf
//
//	script:		
//		
//		name = metadata["name"]
//		basename = name + ".up_matching"
//
//		"""
//		python3 $script $csv $basename
//		"""
//}

