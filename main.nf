#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths

///////////////////////////////////////////////////////////////////////////////
//// METHODS //////////////////////////////////////////////////////////////////

/////////////////////////////////
def addValue(map, key, value) {//
/////////////////////////////////
	def new_map = map.clone()
	new_map.put(key, value)
	return new_map
}

/////////////////////////////
def removeKeys(map, keys) {//
/////////////////////////////
	def new_map = [:]
	map.each{
		if ( ! keys.contains(it.key) )
		{
			new_map.put(it.key, it.value)
		}
	}
	return new_map
}

///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////

include { star } from "./modules/process/align"

include { bam_index as bam_index_tag } from "./modules/process/bam"
include { bam_index as bam_index_select } from "./modules/process/bam"
include { bam_metrics as reads_per_barcode } from "./modules/process/bam"
include { bam_metrics as reads_per_barcode_umi } from "./modules/process/bam"
include { bam_metrics as reads_per_function } from "./modules/process/bam"
include { bam_metrics_hmem as reads_up_matching } from "./modules/process/bam"
include { bam_metrics as reads_barcode_matching } from "./modules/process/bam"
include { bam_metrics as count_duplicates } from "./modules/process/bam"
include { bam_metrics as count_mappings } from "./modules/process/bam"
include { bam_metrics as count_resolved } from "./modules/process/bam"
include { bam_metrics as count_select } from "./modules/process/bam"

include { mark_duplicates } from "./modules/process/tagging"
include { gene } from "./modules/process/tagging"
include { bam_tag as tag_bam } from "./modules/process/tagging"
include { bam_tag_hmem as primary } from "./modules/process/tagging"
include { bam_tag_hmem as multimap } from "./modules/process/tagging"
include { bam_tag_hmem as select } from "./modules/process/tagging"

include { bcl2fastq } from "./modules/process/demultiplexing"
include { merge_lanes } from "./modules/process/demultiplexing"

include { get_barcodes } from "./modules/process/integration"
include { hamming } from "./modules/process/integration"
include { matcher } from "./modules/process/integration"
include { add_match } from "./modules/process/integration"

include { shuffling } from "./modules/process/pucks"

include { fastqc } from "./modules/process/quality_control"
include { reads_per_barcode_fastq } from "./modules/process/quality_control"

include { extract_barcode } from "./modules/process/up_primer"

include { plot_1_arg as plot_balance_barcode } from "./modules/process/plot"
include { plot_1_arg as plot_balance_umi } from "./modules/process/plot"
include { plot_1_arg as plot_barcode_align } from "./modules/process/plot"
include { plot_1_arg as plot_barcode_extraction } from "./modules/process/plot"
include { plot_1_arg as plot_barcode_matching } from "./modules/process/plot"
include { plot_1_arg as plot_duplicates } from "./modules/process/plot"
include { plot_1_arg as plot_histo_errors } from "./modules/process/plot"
include { plot_1_arg as plot_histo_function } from "./modules/process/plot"
include { plot_1_arg as plot_histo_genes } from "./modules/process/plot"
include { plot_1_arg as plot_histo_umis } from "./modules/process/plot"
include { plot_1_arg as plot_mappings } from "./modules/process/plot"
include { plot_1_arg as plot_reads_fraction } from "./modules/process/plot"
include { plot_1_arg as plot_resolved } from "./modules/process/plot"
include { plot_1_arg as plot_select } from "./modules/process/plot"
include { plot_1_arg as plot_umis_per_barcode } from "./modules/process/plot"
include { plot_1_arg as plot_up_align } from "./modules/process/plot"
include { plot_1_arg as plot_up_matching } from "./modules/process/plot"
include { plot_2_args as plot_histo_hamming } from "./modules/process/plot"
include { plot_2_args as plot_spatial_umis } from "./modules/process/plot"

include { merge_plots } from "./modules/process/export"
include { rename_coords } from "./modules/process/export"
include { dge } from "./modules/process/export"

///////////////////////////////////////////////////////////////////////////////
//// SCRIPTS //////////////////////////////////////////////////////////////////

extract_barcode_script = Channel.fromPath("bin/extract_barcode")
tag_bam_script = Channel.fromPath("bin/tag_bam")
shuffling_script = Channel.fromPath("bin/shuffling.py")
matcher_script = Channel.fromPath("bin/matcher.py")
add_match_script = Channel.fromPath("bin/add_match")
gene_script = Channel.fromPath("bin/gene")
primary_script = Channel.fromPath("bin/primary")
multimap_script = Channel.fromPath("bin/multimap")
select_script = Channel.fromPath("bin/select")
count_script = Channel.fromPath("bin/count")

reads_per_barcode_fastq_script = Channel.fromPath("bin/fastq/reads_per_barcode.py")

reads_per_barcode_script = Channel.fromPath("bin/bam/reads_per_barcode.py")
reads_up_matching_script = Channel.fromPath("bin/bam/reads_up_matching.py")
reads_barcode_matching_script = Channel.fromPath("bin/bam/reads_barcode_matching.py")
reads_per_barcode_umi_script = Channel.fromPath("bin/bam/reads_per_barcode_umi.py")
count_duplicates_script = Channel.fromPath("bin/bam/count_duplicates.py")
count_mappings_script = Channel.fromPath("bin/bam/count_mappings.py")
count_resolved_script = Channel.fromPath("bin/bam/count_resolved.py")
count_select_script = Channel.fromPath("bin/bam/count_select.py")
reads_per_function_script = Channel.fromPath("bin/bam/reads_per_function.py")

hamming_script =
	Channel
		.fromPath("bin/hamming/hamming")
		.concat(
			Channel
				.fromPath("bin/hamming/cl")
		)
		.collect()

// plots
plot_balance_barcode_script = Channel.fromPath("bin/plot/balance_barcode.py")
plot_balance_umi_script  = Channel.fromPath("bin/plot/balance_umi.py")
plot_barcode_align_script = Channel.fromPath("bin/plot/barcode_align.py")
plot_barcode_extraction_script  = Channel.fromPath("bin/plot/barcode_extraction.py")
plot_barcode_matching_script = Channel.fromPath("bin/plot/barcode_matching.py")
plot_duplicates_script = Channel.fromPath("bin/plot/duplicates.py")
plot_histo_errors_script = Channel.fromPath("bin/plot/histo_errors.py")
plot_histo_function_script = Channel.fromPath("bin/plot/histo_function.py")
plot_histo_genes_script = Channel.fromPath("bin/plot/histo_genes.py")
plot_histo_hamming_script = Channel.fromPath("bin/plot/histo_hamming.py")
plot_histo_umis_script = Channel.fromPath("bin/plot/histo_umis.py")
plot_mappings_script = Channel.fromPath("bin/plot/mappings.py")
plot_reads_fraction_script = Channel.fromPath("bin/plot/reads_fraction.py")
plot_resolved_script = Channel.fromPath("bin/plot/resolved.py")
plot_select_script = Channel.fromPath("bin/plot/select.py")
plot_spatial_umis_script  = Channel.fromPath("bin/plot/spatial_umi.py")
plot_umis_per_barcode_script = Channel.fromPath("bin/plot/umis_per_barcode.py")
plot_up_align_script = Channel.fromPath("bin/plot/up_align.py")
plot_up_matching_script = Channel.fromPath("bin/plot/up_matching.py")

///////////////////////////////////////////////////////////////////////////////
//// DESIGN ///////////////////////////////////////////////////////////////////

Channel
	.fromPath(params.design)
	.splitCsv(header: true)
	.set{ FASTQ }

///////////////////////////////////////////////////////////////////////////////
//// PUCKS ////////////////////////////////////////////////////////////////////

Channel
	.fromPath("pucks/*")
	.map{[ it.getFileName().toString().replaceAll("(\\S+)\\.csv", "\$1"), it ]}
	.set{ PUCKS }

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////

workflow {

	///////////////////////////////////////////////////////////////////////////
	// MERGE

	FASTQ
		.map{[
			it["name"],
			removeKeys(it, ["fastq_1", "fastq_2"]),
			it["fastq_1"],
			it["fastq_2"]
		]}
		.groupTuple()
		.map{[
			it[1][0],
			it[2].sort{ new File(it).getName() },
			it[3].sort{ new File(it).getName() }
		]}
		.set{ TO_MERGE }
	
	merge_lanes(TO_MERGE)

	///////////////////////////////////////////////////////////////////////////

	FASTQ
		.map{[
			removeKeys(it, ["fastq_1", "fastq_2"]),
			it["fastq_1"],
			it["fastq_2"]
		]}
		.set{ TO_FASTQC }

	fastqc(TO_FASTQC)

	////////////////////////////////////////////////////////////////////////////
	// PUCK BARCODES

	shuffling( PUCKS.combine(shuffling_script) )

	shuffling
		.out
		.ordered
		.concat( shuffling.out.shuffled )
		.set{ PUCK_BARCODES }

	///////////////////////////////////////////////////////////////////////////
	// SEQUENCING BARCODES

	extract_barcode( merge_lanes.out.combine(extract_barcode_script) )

	///////////////////////////////////////////////////////////////////////////
	// ALIGNMENT

	star(extract_barcode.out.fastq)

	///////////////////////////////////////////////////////////////////////////
	// BAM FILE TAGGING

	mark_duplicates(star.out.bam)

	tag_bam(
		mark_duplicates
			.out
			.bam
			.combine( Channel.from("tagged") )
			.combine(tag_bam_script)
	)

	bam_index_tag(tag_bam.out)

	reads_per_barcode(
		bam_index_tag
			.out
			.combine( Channel.from("reads_per_barcode") )
			.combine(reads_per_barcode_script)
	)
	get_barcodes( reads_per_barcode.out )

	///////////////////////////////////////////////////////////////////////////
	// MATCHING

	get_barcodes
		.out
		.combine( PUCK_BARCODES )
		.filter{ it[0]["puck"] == it[2] }
		.map{ [ addValue(it[0], "barcodes", it[3]) , it[1] , it[4] ] }
		.set{ TO_HAMMING }

	hamming( TO_HAMMING.combine(hamming_script) )

	hamming
		.out
		.combine(reads_per_barcode.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ it[0] , it[1] , it[3] ] }
		.combine( PUCKS )
		.filter{ it[0]["puck"] == it[3] }
		.map{ [ *it[0..2] , it[4] ] }
		.set{ TO_MATCHING }
	
	matcher( TO_MATCHING.combine(matcher_script) )

	matcher
		.out
		.mapping
		.filter{ it[0]["barcodes"] == "ordered" }
		.combine(tag_bam.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ *it[0..1] , it[3] ] }
		.set{ TO_ADD_MATCH }

	add_match( TO_ADD_MATCH.combine(add_match_script) )

	///////////////////////////////////////////////////////////////////////////
	// GENE COUNT

	gene( add_match.out.combine(gene_script) )
	primary(
		gene
			.out
			.combine( Channel.from("primary") )
			.combine(primary_script)
	)
	multimap(
		primary
			.out
			.combine( Channel.from("multimap") )
			.combine(multimap_script)
	)
	select(
		multimap
			.out
			.combine( Channel.from("select") )
			.combine(select_script)
	)

	dge( select.out.combine(count_script) )

	bam_index_select(select.out)

	///////////////////////////////////////////////////////////////////////////
	// METRICS

	reads_per_barcode_fastq(
		extract_barcode.out.fastq.combine(reads_per_barcode_fastq_script)
	)

	reads_up_matching(
		bam_index_select
		.out
		.combine( Channel.from("reads_up_matching") )
		.combine(reads_up_matching_script)
	)
	reads_barcode_matching(
		bam_index_select
		.out
		.combine( Channel.from("reads_barcode_matching") )
		.combine(reads_barcode_matching_script)
	)
	reads_per_barcode_umi(
		bam_index_select
		.out
		.combine( Channel.from("reads_per_barcode_umi") )
		.combine(reads_per_barcode_umi_script)
	)
	count_duplicates(
		bam_index_select
		.out
		.combine( Channel.from("count_duplicates") )
		.combine(count_duplicates_script)
	)
	count_mappings(
		bam_index_select
		.out
		.combine( Channel.from("count_mappings") )
		.combine(count_mappings_script)
	)
	count_resolved(
		bam_index_select
		.out
		.combine( Channel.from("count_resolved") )
		.combine(count_resolved_script)
	)
	count_select(
		bam_index_select
		.out
		.combine( Channel.from("count_select") )
		.combine(count_select_script)
	)
	reads_per_function(
		bam_index_select
		.out
		.combine( Channel.from("reads_per_function") )
		.combine(reads_per_function_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// PLOTS

	plot_balance_barcode(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("balance_barcode") )
			.combine(plot_balance_barcode_script)
	)
	plot_balance_umi(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("balance_umi") )
			.combine(plot_balance_umi_script)
	)
	plot_barcode_align(
		reads_barcode_matching
			.out
			.combine( Channel.from("barcode_align") )
			.combine(plot_barcode_align_script)
	)
	plot_barcode_extraction(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from("barcode_extraction") )
			.combine(plot_barcode_extraction_script)
	)
	plot_duplicates(
		count_duplicates
			.out
			.combine( Channel.from("duplicates") )
			.combine(plot_duplicates_script)
	)
	plot_mappings(
		count_mappings
			.out
			.combine( Channel.from("mappings") )
			.combine(plot_mappings_script)
	)
	plot_histo_function(
		reads_per_function
			.out
			.combine( Channel.from("histo_function") )
			.combine(plot_histo_function_script)
	)
	plot_reads_fraction(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("reads_fraction") )
			.combine(plot_reads_fraction_script)
	)
	plot_resolved(
		count_resolved
			.out
			.combine( Channel.from("resolved") )
			.combine(plot_resolved_script)
	)
	plot_select(
		count_select
			.out
			.combine( Channel.from("selected") )
			.combine(plot_select_script)
	)
	plot_up_matching(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from("up_matching") )
			.combine(plot_up_matching_script)
	)
	plot_up_align(
		reads_up_matching
			.out
			.combine( Channel.from("up_align") )
			.combine(plot_up_align_script)
	)
	plot_barcode_matching(
		matcher
			.out
			.metrics
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("barcode_matching") )
			.combine(plot_barcode_matching_script)
	)
	plot_histo_errors(
		matcher
			.out
			.mapping
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("histo_errors") )
			.combine(plot_histo_errors_script)
	)

	// with dge/mtx
	plot_histo_umis(
		dge
			.out
			.combine( Channel.from("histo_umis") )
			.combine(plot_histo_umis_script)
	)
	plot_histo_genes(
		dge
			.out
			.combine( Channel.from("histo_genes") )
			.combine(plot_histo_genes_script)
	)
	plot_umis_per_barcode(
		dge
			.out
			.combine( Channel.from("umis_per_barcode") )
			.combine(plot_umis_per_barcode_script)
	)

	////////////////
	// histo hamming
	plot_histo_hamming(
		hamming
			.out
			.map{ [ it[0]["name"] , *it ] }
			.groupTuple()
			.map{
				it[1][0]["barcodes"] == "ordered" ? [it[1][0], *it[2]] : [it[1][1], *it[2].reverse()]
			}
			.combine( Channel.from("histo_hamming") )
			.combine(plot_histo_hamming_script)
	)
	////////////////

	///////////////
	// spatial umis
	plot_spatial_umis(
		dge
			.out
			.combine(matcher.out.coords)
			.filter{ it[2]["barcodes"] == "ordered" }
			.filter{ it[0]["name"] == it[2]["name"] }
			.map{ [ * it[0..1] , it[3] ] }
			.combine( Channel.from("spatial_umis") )
			.combine(plot_spatial_umis_script)
	)
	////////////////

	///////////////////////////////////////////////////////////////////////////
	// REPORT
	
	plot_barcode_extraction
		.out
		.pdf
		.concat(
			plot_up_matching.out.pdf,
			plot_up_align.out.pdf,
			plot_barcode_align.out.pdf,
			plot_duplicates.out.pdf,
			plot_mappings.out.pdf,
			plot_resolved.out.pdf,
			plot_select.out.pdf,
			plot_histo_function.out.pdf,
			plot_balance_barcode.out.pdf,
			plot_balance_umi.out.pdf,
			plot_reads_fraction.out.pdf,
			plot_histo_umis.out.pdf,
			plot_histo_genes.out.pdf,
			plot_barcode_matching.out.pdf,
			plot_histo_errors.out.pdf,
			plot_histo_hamming.out.pdf,
			plot_umis_per_barcode.out.pdf,
			plot_spatial_umis.out.pdf
		)
		.map{ [ it[0]["name"] , it[1] ] }
		.groupTuple()
		.set{ TO_MERGE_PLOTS }
	
	merge_plots(TO_MERGE_PLOTS)
	rename_coords( matcher.out.coords.filter{ it[0]["barcodes"] == "ordered"} )
}

