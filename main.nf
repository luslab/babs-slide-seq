#!/usr/bin/env nextflow

import java.nio.file.Paths
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl=2

///////////////////////////////////////////////////////////////////////////////
//// FUNCTONS /////////////////////////////////////////////////////////////////

//
// NOURDINE CUSTOM FUNCTIONS
//

def absPath(path) {
	def f = new File(path)

	if ( ! f.isAbsolute() ) {
		return Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		return path
	}
}

def addValue(map, key, value) {
	def new_map = map.clone()
	new_map.put(key, value)
	return new_map
}

def parseSeries(metadata, path) {
	def csv = parseCsv( new File(path).text )
	def channels = []
	for ( row in csv ) {
		channels.add( metadata.clone() << row.toMap() )
	}
	return channels
}

def dropKeys(map, keys) {
	// because Map.dropWhile doesn't with the current Java version of Nextflow,
	// apparently requires Java9
	def new_map = [:]
	map.each{ k, v ->
		if ( ! keys.contains(k) ) {
			new_map[k] = v
		}
	}
	return new_map
}

def removeKeys(map, keys) {
	def new_map = [:]

	map.each{
		if ( ! keys.contains(it.key) )
		{
			new_map.put(it.key, it.value)
		}
	}

	return new_map
}

def getMinLength(structure) {
	return structure.split("[A-Z]").collect{it as int }.sum()
}

def getPuckName(puck) {
	def f = new File(puck)
	return f.getName().toString().replaceAll('\\.csv$', '') // single quotes!
}

///////////////////////////////////////////////////////////////////////////////
//// PARAMETER CHECKING ///////////////////////////////////////////////////////

params.fasta        = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.star_index   = WorkflowMain.getGenomeAttribute(params, 'star')
params.gtf          = WorkflowMain.getGenomeAttribute(params, 'gtf')
// params.gene_bed  = WorkflowMain.getGenomeAttribute(params, 'bed12')
// params.blacklist = WorkflowMain.getGenomeAttribute(params, 'blacklist')

checkPathParamList = [
	params.design,
	params.fasta,
    params.star_index,
    params.gtf,
	// params.blacklist
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.design) { ch_input = file(params.design) } else { exit 1, "Design samplesheet not specified!" }

// ch_blacklist = Channel.empty()
// if (params.blacklist) {
//     ch_blacklist = file(params.blacklist)
// }
// else {
//     ch_blacklist = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
//     WorkflowCutandrun.blacklistWarn(log)
// }

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// // Stage dummy file to be used as an optional input where required
// ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)


///////////////////////////////////////////////////////////////////////////////
//// INCLUDES ////////////////////////////////////////////////////////////////


//////////
// samples

include { BCL2FASTQ   } from "./modules/process/demultiplexing"
include { MERGE_LANES } from "./modules/process/demultiplexing"
include { FASTQC      } from "./modules/process/quality_control"
//////////

////////
// pucks
include { SHUFFLING } from "./modules/process/pucks"
////////

/////////////////////
// barcode extraction

include { EXTRACT_BARCODE } from "./modules/process/extract_barcode"
include { PLOT as PLOT_UP_MATCHING } from "./modules/process/plot"
include { PLOT as PLOT_BARCODE_EXTRACTION } from "./modules/process/plot"
/////////////////////

// ///////////////////////////
// // alignment and duplicates
// include { star } from "./modules/process/align"
// include { mark_duplicates } from "./modules/process/tagging"
// ///////////////////////////

// /////////////////
// // slide seq tags
// include { bam_tag as tag_bam } from "./modules/process/tagging"
// tag_bam_script = Channel.fromPath("bin/tag_bam")

// include { bam_metrics_hmem as reads_up_matching } from "./modules/process/bam"
// reads_up_matching_script = Channel.fromPath("bin/bam/reads_up_matching.py")

// include { plot_1_arg as plot_up_align } from "./modules/process/plot"
// plot_up_align_script = Channel.fromPath("bin/plot/up_align.py")

// include { bam_filter as bam_filter_up_matched } from "./modules/process/bam"
// /////////////////

// /////////////////////////////
// // umis per barcode threshold

// include { umis_per_barcode } from "./modules/process/tagging"
// umis_per_barcode_script = Channel.fromPath("bin/umis_per_barcode")

// include { bam_metrics as reads_umis_per_barcode } from "./modules/process/bam"
// reads_umis_per_barcode_script = Channel.fromPath("bin/bam/reads_umis_per_barcode.py")

// include { bam_metrics as reads_umi_threshold } from "./modules/process/bam"
// reads_umi_threshold_script = Channel.fromPath("bin/bam/reads_umi_threshold.py")

// include { bam_filter as bam_filter_umi_threshold } from "./modules/process/bam"

// include { bam_metrics as reads_barcode_matching } from "./modules/process/bam"
// reads_barcode_matching_script = Channel.fromPath("bin/bam/reads_barcode_matching.py")

// include { plot_1_val as plot_umi_threshold } from "./modules/process/plot"
// plot_umi_threshold_script = Channel.fromPath("bin/plot/umi_threshold.py")
// /////////////////////////////

// ///////////////////
// // hamming distance

// include { hamming } from "./modules/process/integration"
// hamming_script =
// 	Channel
// 		.fromPath("bin/hamming/hamming")
// 		.concat(
// 			Channel
// 				.fromPath("bin/hamming/cl")
// 		)
// 		.collect()

// include { plot_2_args as plot_histo_hamming } from "./modules/process/plot"
// plot_histo_hamming_script = Channel.fromPath("bin/plot/histo_hamming.py")
// ///////////////////

// ///////////////////
// // barcode matching

// include { get_barcodes } from "./modules/process/integration"

// include { matcher } from "./modules/process/integration"
// matcher_script = Channel.fromPath("bin/matcher.py")

// include { add_match } from "./modules/process/integration"
// add_match_script = Channel.fromPath("bin/add_match")

// include { plot_1_arg as plot_barcode_matching } from "./modules/process/plot"
// plot_barcode_matching_script = Channel.fromPath("bin/plot/barcode_matching.py")

// include { plot_1_arg as plot_barcode_align } from "./modules/process/plot"
// plot_barcode_align_script = Channel.fromPath("bin/plot/barcode_align.py")

// include { plot_1_arg as plot_histo_errors } from "./modules/process/plot"
// plot_histo_errors_script = Channel.fromPath("bin/plot/histo_errors.py")

// include { bam_filter as bam_filter_barcode_matched } from "./modules/process/bam"
// ///////////////////

// ///////////////
// // gene tagging

// include { htseq } from "./modules/process/tagging"

// include { bam_metrics as count_gene_tags } from "./modules/process/bam"
// count_gene_tags_script = Channel.fromPath("bin/bam/count_gene_tags.py")

// include { bam_metrics as count_reads_per_umi } from "./modules/process/bam"
// count_reads_per_umi_script = Channel.fromPath("bin/bam/reads_per_umi.py")

// include { bam_metrics as count_reads_per_umi_gene } from "./modules/process/bam"
// count_reads_per_umi_gene_script = Channel.fromPath("bin/bam/reads_per_umi_gene.py")

// include { plot_1_arg as plot_gene_tags } from "./modules/process/plot"
// plot_gene_tags_script = Channel.fromPath("bin/plot/gene_tags.py")

// include { bam_filter as bam_filter_gene_tags } from "./modules/process/bam"

// ///////////////

// /////////////////////
// // umis multi mapping

// include { select } from "./modules/process/tagging"
// select_script = Channel.fromPath("bin/select")

// include { duplicates } from "./modules/process/quality_control"
// duplicates_script = Channel.fromPath("bin/duplicates.py")

// include { bam_metrics as count_select } from "./modules/process/bam"
// count_select_script = Channel.fromPath("bin/bam/count_select.py")

// include { plot_1_arg as plot_select } from "./modules/process/plot"
// plot_select_script = Channel.fromPath("bin/plot/select.py")

// include { bam_filter as bam_filter_multimapped_umis } from "./modules/process/bam"
// /////////////////////

// ////////////
// // sequences

// include { bam_metrics as reads_per_barcode_umi } from "./modules/process/bam"
// reads_per_barcode_umi_script = Channel.fromPath("bin/bam/reads_per_barcode_umi.py")

// include { plot_1_arg as plot_balance_barcode } from "./modules/process/plot"
// plot_balance_barcode_script = Channel.fromPath("bin/plot/balance_barcode.py")

// include { plot_1_arg as plot_balance_umi } from "./modules/process/plot"
// plot_balance_umi_script  = Channel.fromPath("bin/plot/balance_umi.py")

// include { plot_1_arg as plot_reads_fraction } from "./modules/process/plot"
// plot_reads_fraction_script = Channel.fromPath("bin/plot/reads_fraction.py")
// ////////////

// //////
// // dge

// include { dge } from "./modules/process/export"
// count_script = Channel.fromPath("bin/count")

// plot_histo_genes_script = Channel.fromPath("bin/plot/histo_genes.py")
// plot_histo_umis_script = Channel.fromPath("bin/plot/histo_umis.py")
// plot_spatial_umis_script  = Channel.fromPath("bin/plot/spatial_umi.py")
// plot_umis_per_barcode_script = Channel.fromPath("bin/plot/umis_per_barcode.py")

// include { plot_1_arg as plot_histo_genes } from "./modules/process/plot"
// include { plot_1_arg as plot_histo_umis } from "./modules/process/plot"
// include { plot_1_arg as plot_umis_per_barcode } from "./modules/process/plot"
// include { plot_2_args as plot_spatial_umis } from "./modules/process/plot"
// //////

// /////////
// // output 

// include { merge_plots } from "./modules/process/export"
// include { rename_coords } from "./modules/process/export"
// /////////


///////////////////////////////////////////////////////////////////////////////
//// DESIGN ///////////////////////////////////////////////////////////////////

Channel
	.fromPath(params.design)
	.splitCsv(header: true)
	.map{ addValue(it, "min_length", getMinLength(it["read_structure"])) }
	.map{ addValue(it, "puck_path", new File(it["puck"]).getAbsolutePath()) }
	.map{ addValue(it, "puck", getPuckName(it["puck_path"])) }
	.map{ addValue(it, "fastq_1", new File(it["fastq_1"]).getAbsolutePath()) }
	.map{ addValue(it, "fastq_2", new File(it["fastq_2"]).getAbsolutePath()) }
	.set{ ch_fastq }
//ch_fastq | view

///////////////////////////////////////////////////////////////////////////////
//// PUCKS ////////////////////////////////////////////////////////////////////

ch_fastq
	.map{ [ it["puck"] , it["puck_path"] ] }
	.unique()
	.set{ ch_pucks }
//ch_pucks | view

///////////////////////////////////////////////////////////////////////////////
//// MAIN WORKFLOW ////////////////////////////////////////////////////////////

workflow {
	///////////////////////////////////////////////////////////////////////////
	// MERGE

	ch_fastq
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
		.set{ ch_fastq_merge }
	
	MERGE_LANES( ch_fastq_merge )

	///////////////////////////////////////////////////////////////////////////

	ch_fastq
		.map{[
			removeKeys(it, ["fastq_1", "fastq_2"]),
			it["fastq_1"],
			it["fastq_2"]
		]}
		.set{ ch_fastqc }

	FASTQC( ch_fastqc )

	///////////////////////////////////////////////////////////////////////////
	// PUCK BARCODES

	SHUFFLING( ch_pucks )

	SHUFFLING
		.out
		.ordered
		.concat( SHUFFLING.out.shuffled )
		.set{ ch_puck_barcodes }
	//ch_puck_barcodes | view

	///////////////////////////////////////////////////////////////////////////
	// SEQUENCING BARCODES

	EXTRACT_BARCODE( MERGE_LANES.out )
	//EXTRACT_BARCODE.out.metrics | view

	PLOT_UP_MATCHING( EXTRACT_BARCODE.out.metrics.map{ [ it[0] , it[1], '' ] } )
	PLOT_BARCODE_EXTRACTION( EXTRACT_BARCODE.out.metrics.map{ [ it[0] , it[1], it[0]["min_length"] ] } )

	// ///////////////////////////////////////////////////////////////////////////
	// // ALIGNMENT

	// star(extract_barcode.out.fastq)

	// ///////////////////////////////////////////////////////////////////////////
	// // DUPLICATES

	// mark_duplicates(star.out.bam)

	// ///////////////////////////////////////////////////////////////////////////
	// // ADD SLIDE-SEQ, ALIGNMENT AND DUPLICATES TAGS

	// // slide-seq tags, alignment tag and duplicate tag
	// tag_bam(
	// 	mark_duplicates
	// 		.out
	// 		.bam
	// 		.combine( Channel.from("tagged") )
	// 		.combine(tag_bam_script)
	// )

	// // 3 columns: Matched, Mapped, Reads
	// reads_up_matching(
	// 	tag_bam
	// 		.out
	// 		.combine( Channel.from("reads_up_matching") )
	// 		.combine(reads_up_matching_script)
	// )
	// plot_up_align(
	// 	reads_up_matching
	// 		.out
	// 		.combine( Channel.from("up_align") )
	// 		.combine(plot_up_align_script)
	// )

	// bam_filter_up_matched(
	// 	tag_bam
	// 		.out
	// 		.combine( Channel.from("up_matched") )
	// 		.combine( Channel.from("[us]==\"MATCHED\" && [as]==\"MAPPED\"") )
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // UMIS PER BARCODE THRESHOLD

	// umis_per_barcode(
	// 	bam_filter_up_matched
	// 		.out
	// 		.combine(umis_per_barcode_script)
	// )

	// reads_umis_per_barcode(
	// 	umis_per_barcode
	// 		.out
	// 		.combine( Channel.from("reads_umis_per_barcode") )
	// 		.combine(reads_umis_per_barcode_script)
	// )

	// reads_umi_threshold(
	// 	umis_per_barcode
	// 		.out
	// 		.combine( Channel.from("reads_umi_threshold") )
	// 		.combine(reads_umi_threshold_script)
	// )

	// plot_umi_threshold(
	// 	reads_umi_threshold
	// 		.out
	// 		.combine( Channel.from(params.umis_threshold) )
	// 		.combine( Channel.from("umi_threshold") )
	// 		.combine(plot_umi_threshold_script)
	// )

	// bam_filter_umi_threshold(
	// 	umis_per_barcode
	// 		.out
	// 		.combine( Channel.from("umi_threshold") )
	// 		.combine( Channel.from("[bt]==\"PASS\"") )
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // HAMMING DISTANCE

	// get_barcodes( reads_umis_per_barcode.out )

	// get_barcodes
	// 	.out
	// 	.combine( PUCK_BARCODES )
	// 	.filter{ it[0]["puck"] == it[2] }
	// 	.map{ [ addValue(it[0], "barcodes", it[3]) , it[1] , it[4] ] }
	// 	.set{ TO_HAMMING }

	// hamming( TO_HAMMING.combine(hamming_script) )

	// plot_histo_hamming(
	// 	hamming
	// 		.out
	// 		.map{ [ it[0]["name"] , *it ] }
	// 		.groupTuple()
	// 		.map{
	// 			it[1][0]["barcodes"] == "ordered" ? [it[1][0], *it[2]] : [it[1][1], *it[2].reverse()]
	// 		}
	// 		.combine( Channel.from("histo_hamming") )
	// 		.combine(plot_histo_hamming_script)
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // BARCODE MATCHING

	// hamming
	// 	.out
	// 	.combine(reads_umis_per_barcode.out)
	// 	.filter{ it[0]["name"] == it[2]["name"] }
	// 	.map{ [ it[0] , it[1] , it[3] ] }
	// 	.combine( PUCKS )
	// 	.filter{ it[0]["puck"] == it[3] }
	// 	.map{ [ *it[0..2] , it[4] ] }
	// 	.set{ TO_MATCHING }
	
	// matcher( TO_MATCHING.combine(matcher_script) )

	// matcher
	// 	.out
	// 	.mapping
	// 	.filter{ it[0]["barcodes"] == "ordered" }
	// 	.combine(bam_filter_umi_threshold.out)
	// 	.filter{ it[0]["name"] == it[2]["name"] }
	// 	.map{ [ *it[0..1] , it[3] ] }
	// 	.set{ TO_ADD_MATCH }

	// plot_barcode_matching(
	// 	matcher
	// 		.out
	// 		.metrics
	// 		.filter{ it[0]["barcodes"] == "ordered" }
	// 		.combine( Channel.from("barcode_matching") )
	// 		.combine(plot_barcode_matching_script)
	// )

	// plot_histo_errors(
	// 	matcher
	// 		.out
	// 		.mapping
	// 		.filter{ it[0]["barcodes"] == "ordered" }
	// 		.combine( Channel.from("histo_errors") )
	// 		.combine(plot_histo_errors_script)
	// )

	// add_match( TO_ADD_MATCH.combine(add_match_script) )

	// reads_barcode_matching(
	// 	add_match
	// 		.out
	// 		.combine( Channel.from("reads_barcode_matching") )
	// 		.combine(reads_barcode_matching_script)
	// )

	// plot_barcode_align(
	// 	reads_barcode_matching
	// 		.out
	// 		.combine( Channel.from("barcode_align") )
	// 		.combine(plot_barcode_align_script)
	// )

	// bam_filter_barcode_matched(
	// 	add_match
	// 		.out
	// 		.combine( Channel.from("barcode_matched") )
	// 		.combine( Channel.from("[bs]==\"MATCHED\"") )
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // GENE COUNT

	// htseq(bam_filter_barcode_matched.out)

	// count_gene_tags(
	// 	htseq
	// 		.out
	// 		.bam
	// 		.combine( Channel.from("gene_tags") )
	// 		.combine(count_gene_tags_script)
	// )

	// plot_gene_tags(
	// 	count_gene_tags
	// 		.out
	// 		.combine( Channel.from("gene_tags") )
	// 		.combine(plot_gene_tags_script)
	// )

	// bam_filter_gene_tags(
	// 	htseq
	// 		.out
	// 		.bam
	// 		.combine( Channel.from("gene_tags") )
	// 		.combine( Channel.from("[XF]!~\"^__.+\"") )
	// )

	// count_reads_per_umi(
	// 	bam_filter_gene_tags
	// 		.out
	// 		.combine( Channel.from("reads_per_umi") )
	// 		.combine(count_reads_per_umi_script)
	// )

	// count_reads_per_umi_gene(
	// 	bam_filter_gene_tags
	// 		.out
	// 		.combine( Channel.from("reads_per_umi_gene") )
	// 		.combine(count_reads_per_umi_gene_script)
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // UMIS MAPPINGS

	// select(
	// 	bam_filter_gene_tags
	// 		.out
	// 		.map{ it[0..1] }
	// 		.combine( Channel.from("select") )
	// 		.combine(select_script)
	// )

	// select
	// 	.out
	// 	.unique_reads
	// 	.concat(select.out.resolved_reads)
	// 	.concat(select.out.unresolved_reads)
	// 	.combine(merge_lanes.out)
	// 	.filter{ it[0]["name"] == it[3]["name"] }
	// 	.map{ [ addValue(it[0], "status", it[1]) , it[2] , *it[4..5] ] }
	// 	.set{ TO_DUPLICATES }
	// duplicates(TO_DUPLICATES.combine(duplicates_script))

	// count_select(
	// 	select
	// 		.out
	// 		.bam
	// 		.combine( Channel.from("count_select") )
	// 		.combine(count_select_script)
	// )
	// plot_select(
	// 	count_select
	// 		.out
	// 		.combine( Channel.from("selected") )
	// 		.combine(plot_select_script)
	// )

	// bam_filter_multimapped_umis(
	// 	select
	// 		.out
	// 		.bam
	// 		.combine( Channel.from("multimap_umis") )
	// 		.combine( Channel.from("[cs]==\"UNIQUE\" || [cs]==\"INCLUDED\"") )
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // SEQUENCES

	// reads_per_barcode_umi(
	// 	bam_filter_multimapped_umis
	// 		.out
	// 		.combine( Channel.from("reads_per_barcode_umi") )
	// 		.combine(reads_per_barcode_umi_script)
	// )

	// plot_balance_barcode(
	// 	reads_per_barcode_umi
	// 		.out
	// 		.combine( Channel.from("balance_barcode") )
	// 		.combine(plot_balance_barcode_script)
	// )

	// plot_balance_umi(
	// 	reads_per_barcode_umi
	// 		.out
	// 		.combine( Channel.from("balance_umi") )
	// 		.combine(plot_balance_umi_script)
	// )

	// plot_reads_fraction(
	// 	reads_per_barcode_umi
	// 		.out
	// 		.combine( Channel.from("reads_fraction") )
	// 		.combine(plot_reads_fraction_script)
	// )

	// ///////////////////////////////////////////////////////////////////////////
	// // EXPRESSION MATRIX

	// dge( bam_filter_multimapped_umis.out.map{it[0..1]}.combine(count_script) )

	// plot_histo_umis(
	// 	dge
	// 		.out
	// 		.combine( Channel.from("histo_umis") )
	// 		.combine(plot_histo_umis_script)
	// )

	// plot_histo_genes(
	// 	dge
	// 		.out
	// 		.combine( Channel.from("histo_genes") )
	// 		.combine(plot_histo_genes_script)
	// )

	// plot_umis_per_barcode(
	// 	dge
	// 		.out
	// 		.combine( Channel.from("umis_per_barcode") )
	// 		.combine(plot_umis_per_barcode_script)
	// )

	// ///////////////
	// // spatial umis
	// plot_spatial_umis(
	// 	dge
	// 		.out
	// 		.combine(matcher.out.coords)
	// 		.filter{ it[2]["barcodes"] == "ordered" }
	// 		.filter{ it[0]["name"] == it[2]["name"] }
	// 		.map{ [ * it[0..1] , it[3] ] }
	// 		.combine( Channel.from("spatial_umis") )
	// 		.combine(plot_spatial_umis_script)
	// )
	// ////////////////

	// ///////////////////////////////////////////////////////////////////////////
	// // OUTPUT
	
	// plot_barcode_extraction
	// 	.out
	// 	.pdf
	// 	.concat(
	// 		plot_up_matching.out.pdf,
	// 		plot_up_align.out.pdf,
	// 		plot_umi_threshold.out.pdf,
	// 		plot_barcode_align.out.pdf,
	// 		plot_gene_tags.out.pdf,
	// 		plot_select.out.pdf,
	// 		plot_balance_barcode.out.pdf,
	// 		plot_balance_umi.out.pdf,
	// 		plot_reads_fraction.out.pdf,
	// 		plot_histo_umis.out.pdf,
	// 		plot_histo_genes.out.pdf,
	// 		plot_barcode_matching.out.pdf,
	// 		plot_histo_errors.out.pdf,
	// 		plot_histo_hamming.out.pdf,
	// 		plot_umis_per_barcode.out.pdf,
	// 		plot_spatial_umis.out.pdf
	// 	)
	// 	.map{ [ it[0]["name"] , it[1] ] }
	// 	.groupTuple()
	// 	.set{ TO_MERGE_PLOTS }
	
	// merge_plots(TO_MERGE_PLOTS)

	// rename_coords( matcher.out.coords.filter{ it[0]["barcodes"] == "ordered"} )
}

