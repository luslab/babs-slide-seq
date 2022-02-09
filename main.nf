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


//////////
// samples

include { bcl2fastq } from "./modules/process/demultiplexing"
include { merge_lanes } from "./modules/process/demultiplexing"
include { fastqc } from "./modules/process/quality_control"
//////////

////////
// pucks

include { shuffling } from "./modules/process/pucks"
shuffling_script = Channel.fromPath("bin/shuffling.py")
////////

/////////////////////
// barcode extraction

include { extract_barcode } from "./modules/process/up_primer"
extract_barcode_script = Channel.fromPath("bin/extract_barcode")

include { plot_1_arg as plot_up_matching } from "./modules/process/plot"
plot_up_matching_script = Channel.fromPath("bin/plot/up_matching.py")
include { plot_1_arg_1_val as plot_barcode_extraction } from "./modules/process/plot"
plot_barcode_extraction_script  = Channel.fromPath("bin/plot/barcode_extraction.py")
/////////////////////

///////////////////////////
// alignment and duplicates
include { star } from "./modules/process/align"
include { mark_duplicates } from "./modules/process/tagging"
///////////////////////////

/////////////////
// slide seq tags
include { bam_tag as tag_bam } from "./modules/process/tagging"
tag_bam_script = Channel.fromPath("bin/tag_bam")

include { bam_metrics_hmem as reads_up_matching } from "./modules/process/bam"
reads_up_matching_script = Channel.fromPath("bin/bam/reads_up_matching.py")

include { plot_1_arg as plot_up_align } from "./modules/process/plot"
plot_up_align_script = Channel.fromPath("bin/plot/up_align.py")

include { bam_filter as bam_filter_up_matched } from "./modules/process/bam"
/////////////////

/////////////////////////////
// umis per barcode threshold

include { umis_per_barcode } from "./modules/process/tagging"
umis_per_barcode_script = Channel.fromPath("bin/umis_per_barcode")

include { bam_metrics as reads_umis_per_barcode } from "./modules/process/bam"
reads_umis_per_barcode_script = Channel.fromPath("bin/bam/reads_umis_per_barcode.py")

include { bam_metrics as reads_umi_threshold } from "./modules/process/bam"
reads_umi_threshold_script = Channel.fromPath("bin/bam/reads_umi_threshold.py")

include { bam_filter as bam_filter_umi_threshold } from "./modules/process/bam"

include { bam_metrics as reads_barcode_matching } from "./modules/process/bam"
reads_barcode_matching_script = Channel.fromPath("bin/bam/reads_barcode_matching.py")

include { plot_1_val as plot_umi_threshold } from "./modules/process/plot"
plot_umi_threshold_script = Channel.fromPath("bin/plot/umi_threshold.py")
/////////////////////////////

///////////////////
// hamming distance

include { hamming } from "./modules/process/integration"
hamming_script =
	Channel
		.fromPath("bin/hamming/hamming")
		.concat(
			Channel
				.fromPath("bin/hamming/cl")
		)
		.collect()

include { plot_2_args as plot_histo_hamming } from "./modules/process/plot"
plot_histo_hamming_script = Channel.fromPath("bin/plot/histo_hamming.py")
///////////////////

///////////////////
// barcode matching

include { get_barcodes } from "./modules/process/integration"

include { matcher } from "./modules/process/integration"
matcher_script = Channel.fromPath("bin/matcher.py")

include { add_match } from "./modules/process/integration"
add_match_script = Channel.fromPath("bin/add_match")

include { plot_1_arg as plot_barcode_matching } from "./modules/process/plot"
plot_barcode_matching_script = Channel.fromPath("bin/plot/barcode_matching.py")

include { plot_1_arg as plot_barcode_align } from "./modules/process/plot"
plot_barcode_align_script = Channel.fromPath("bin/plot/barcode_align.py")

include { plot_1_arg as plot_histo_errors } from "./modules/process/plot"
plot_histo_errors_script = Channel.fromPath("bin/plot/histo_errors.py")

include { bam_filter as bam_filter_barcode_matched } from "./modules/process/bam"
///////////////////

///////////////
// gene tagging

include { htseq } from "./modules/process/tagging"

include { bam_metrics as count_gene_tags } from "./modules/process/bam"
count_gene_tags_script = Channel.fromPath("bin/bam/count_gene_tags.py")

include { plot_1_arg as plot_gene_tags } from "./modules/process/plot"
plot_gene_tags_script = Channel.fromPath("bin/plot/gene_tags.py")

include { bam_filter as bam_filter_gene_tags } from "./modules/process/bam"

///////////////

/////////////////////
// umis multi mapping

include { bam_tag_hmem as select } from "./modules/process/tagging"
select_script = Channel.fromPath("bin/select")

include { bam_metrics as count_select } from "./modules/process/bam"
count_select_script = Channel.fromPath("bin/bam/count_select.py")

include { plot_1_arg as plot_select } from "./modules/process/plot"
plot_select_script = Channel.fromPath("bin/plot/select.py")

include { bam_filter as bam_filter_multimapped_umis } from "./modules/process/bam"
/////////////////////

////////////
// sequences

include { bam_metrics as reads_per_barcode_umi } from "./modules/process/bam"
reads_per_barcode_umi_script = Channel.fromPath("bin/bam/reads_per_barcode_umi.py")

include { plot_1_arg as plot_balance_barcode } from "./modules/process/plot"
plot_balance_barcode_script = Channel.fromPath("bin/plot/balance_barcode.py")

include { plot_1_arg as plot_balance_umi } from "./modules/process/plot"
plot_balance_umi_script  = Channel.fromPath("bin/plot/balance_umi.py")

include { plot_1_arg as plot_reads_fraction } from "./modules/process/plot"
plot_reads_fraction_script = Channel.fromPath("bin/plot/reads_fraction.py")
////////////

//////
// dge

include { dge } from "./modules/process/export"
count_script = Channel.fromPath("bin/count")

plot_histo_genes_script = Channel.fromPath("bin/plot/histo_genes.py")
plot_histo_umis_script = Channel.fromPath("bin/plot/histo_umis.py")
plot_spatial_umis_script  = Channel.fromPath("bin/plot/spatial_umi.py")
plot_umis_per_barcode_script = Channel.fromPath("bin/plot/umis_per_barcode.py")

include { plot_1_arg as plot_histo_genes } from "./modules/process/plot"
include { plot_1_arg as plot_histo_umis } from "./modules/process/plot"
include { plot_1_arg as plot_umis_per_barcode } from "./modules/process/plot"
include { plot_2_args as plot_spatial_umis } from "./modules/process/plot"
//////

/////////
// output 

include { merge_plots } from "./modules/process/export"
include { rename_coords } from "./modules/process/export"
/////////


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

	////////////////////////////////////////////////////////////////////////////
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
		.map{ [ addValue(it[0], "n_fastq_read1", it[1].size()) , *it[1..2] ] }
		.map{ [ addValue(it[0], "n_fastq_read2", it[2].size()) , *it[1..2] ] }
		.set{ TO_MERGE }
	
	merge_lanes(TO_MERGE)

	////////////////////////////////////////////////////////////////////////////

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

	////////////////////////////////////////////////////////////////////////////
	// SEQUENCING BARCODES

	extract_barcode( merge_lanes.out.combine(extract_barcode_script) )

	plot_up_matching(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from("up_matching") )
			.combine(plot_up_matching_script)
	)

	plot_barcode_extraction(
		extract_barcode
			.out
			.metrics
			.combine( Channel.from(params.min_length) )
			.combine( Channel.from("barcode_extraction") )
			.combine(plot_barcode_extraction_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// ALIGNMENT

	star(extract_barcode.out.fastq)

	///////////////////////////////////////////////////////////////////////////
	// DUPLICATES

	mark_duplicates(star.out.bam)

	///////////////////////////////////////////////////////////////////////////
	// ADD SLIDE-SEQ, ALIGNMENT AND DUPLICATES TAGS

	// slide-seq tags, alignment tag and duplicate tag
	tag_bam(
		mark_duplicates
			.out
			.bam
			.combine( Channel.from("tagged") )
			.combine(tag_bam_script)
	)

	// 3 columns: Matched, Mapped, Reads
	reads_up_matching(
		tag_bam
			.out
			.combine( Channel.from("reads_up_matching") )
			.combine(reads_up_matching_script)
	)
	plot_up_align(
		reads_up_matching
			.out
			.combine( Channel.from("up_align") )
			.combine(plot_up_align_script)
	)

	bam_filter_up_matched(
		tag_bam
			.out
			.combine( Channel.from("up_matched") )
			.combine( Channel.from("[us]==\"MATCHED\" && [as]==\"MAPPED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// UMIS PER BARCODE THRESHOLD

	umis_per_barcode(
		bam_filter_up_matched
			.out
			.combine(umis_per_barcode_script)
	)

	reads_umis_per_barcode(
		umis_per_barcode
			.out
			.combine( Channel.from("reads_umis_per_barcode") )
			.combine(reads_umis_per_barcode_script)
	)

	reads_umi_threshold(
		umis_per_barcode
			.out
			.combine( Channel.from("reads_umi_threshold") )
			.combine(reads_umi_threshold_script)
	)

	plot_umi_threshold(
		reads_umi_threshold
			.out
			.combine( Channel.from(params.umis_threshold) )
			.combine( Channel.from("umi_threshold") )
			.combine(plot_umi_threshold_script)
	)

	bam_filter_umi_threshold(
		umis_per_barcode
			.out
			.combine( Channel.from("umi_threshold") )
			.combine( Channel.from("[bt]==\"PASS\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// HAMMING DISTANCE

	get_barcodes( reads_umis_per_barcode.out )

	get_barcodes
		.out
		.combine( PUCK_BARCODES )
		.filter{ it[0]["puck"] == it[2] }
		.map{ [ addValue(it[0], "barcodes", it[3]) , it[1] , it[4] ] }
		.set{ TO_HAMMING }

	hamming( TO_HAMMING.combine(hamming_script) )

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

	///////////////////////////////////////////////////////////////////////////
	// BARCODE MATCHING

	hamming
		.out
		.combine(reads_umis_per_barcode.out)
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
		.combine(bam_filter_umi_threshold.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ *it[0..1] , it[3] ] }
		.set{ TO_ADD_MATCH }

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

	add_match( TO_ADD_MATCH.combine(add_match_script) )

	reads_barcode_matching(
		add_match
			.out
			.combine( Channel.from("reads_barcode_matching") )
			.combine(reads_barcode_matching_script)
	)

	plot_barcode_align(
		reads_barcode_matching
			.out
			.combine( Channel.from("barcode_align") )
			.combine(plot_barcode_align_script)
	)

	bam_filter_barcode_matched(
		add_match
			.out
			.combine( Channel.from("barcode_matched") )
			.combine( Channel.from("[bs]==\"MATCHED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// GENE COUNT

	htseq(bam_filter_barcode_matched.out)

	count_gene_tags(
		htseq
			.out
			.bam
			.combine( Channel.from("gene_tags") )
			.combine(count_gene_tags_script)
	)

	plot_gene_tags(
		count_gene_tags
			.out
			.combine( Channel.from("gene_tags") )
			.combine(plot_gene_tags_script)
	)

	bam_filter_gene_tags(
		htseq
			.out
			.bam
			.combine( Channel.from("gene_tags") )
			.combine( Channel.from("[XF]!~\"^__.+\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// UMIS MAPPINGS

	select(
		bam_filter_gene_tags
			.out
			.map{ it[0..1] }
			.combine( Channel.from("select") )
			.combine(select_script)
	)

	count_select(
		select
			.out
			.combine( Channel.from("count_select") )
			.combine(count_select_script)
	)
	plot_select(
		count_select
			.out
			.combine( Channel.from("selected") )
			.combine(plot_select_script)
	)

	bam_filter_multimapped_umis(
		select
			.out
			.combine( Channel.from("multimap_umis") )
			.combine( Channel.from("[cs]==\"UNIQUE\" || [cs]==\"INCLUDED\"") )
	)

	///////////////////////////////////////////////////////////////////////////
	// SEQUENCES

	reads_per_barcode_umi(
		bam_filter_multimapped_umis
			.out
			.combine( Channel.from("reads_per_barcode_umi") )
			.combine(reads_per_barcode_umi_script)
	)

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

	plot_reads_fraction(
		reads_per_barcode_umi
			.out
			.combine( Channel.from("reads_fraction") )
			.combine(plot_reads_fraction_script)
	)

	///////////////////////////////////////////////////////////////////////////
	// EXPRESSION MATRIX

	dge( bam_filter_multimapped_umis.out.map{it[0..1]}.combine(count_script) )

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
	// OUTPUT
	
	plot_barcode_extraction
		.out
		.pdf
		.concat(
			plot_up_matching.out.pdf,
			plot_up_align.out.pdf,
			plot_umi_threshold.out.pdf,
			plot_barcode_align.out.pdf,
			plot_gene_tags.out.pdf,
			plot_select.out.pdf,
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

