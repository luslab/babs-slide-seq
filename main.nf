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
// genome

include { UNTAR as UNTAR_STAR_INDEX } from './modules/nf-core/modules/untar/main'

//////////
// samples

include { BCL2FASTQ   } from "./modules/local/demultiplexing"
include { MERGE_LANES } from "./modules/local/demultiplexing"
include { FASTQC      } from "./modules/local/quality_control"
//////////

////////
// pucks
include { SHUFFLING } from "./modules/local/pucks"
////////

/////////////////////
// barcode extraction

include { EXTRACT_BARCODE                 } from "./modules/local/extract_barcode"
include { PLOT as PLOT_UP_MATCHING        } from "./modules/local/plot"
include { PLOT as PLOT_BARCODE_EXTRACTION } from "./modules/local/plot"
/////////////////////

///////////////////////////
// alignment and duplicates
include { STAR            } from "./modules/local/align"
include { MARK_DUPLICATES } from "./modules/local/tagging"
///////////////////////////

/////////////////
// slide seq tags
include { TAG_BAM                             } from "./modules/local/tagging"
include { BAM_METRICS as READS_UP_MATCHING    } from "./modules/local/bam"
include { PLOT as PLOT_UP_ALIGN               } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_UP_MATCHED } from "./modules/local/bam"
/////////////////

/////////////////////////////
// umis per barcode threshold

include { UMIS_PER_BARCODE                       } from "./modules/local/tagging"
include { BAM_METRICS as READS_UMIS_PER_BARCODE  } from "./modules/local/bam"
include { BAM_METRICS as READS_UMI_THRESHOLD     } from "./modules/local/bam"
include { PLOT as PLOT_UMI_THRESHOLD             } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_UMI_THRESHOLD } from "./modules/local/bam"

// /////////////////////////////

// ///////////////////
// // barcode matching

include { GET_BARCODES } from "./modules/local/integration"

// ///////////////////
// // hamming distance

include { HAMMING } from "./modules/local/integration"
include { PLOT_HAMMING_HISTO } from "./modules/local/plot"

// ///////////////////

include { MATCHER } from "./modules/local/integration"
include { PLOT as PLOT_BARCODE_MATCHING } from "./modules/local/plot"
include { PLOT as PLOT_HISTO_ERROR } from "./modules/local/plot"
include { ADD_MATCH } from "./modules/local/integration"
include { BAM_METRICS as READS_BARCODE_MATCHING } from "./modules/local/bam"
include { PLOT as PLOT_BARCODE_ALIGN } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_BARCODE_MATCHED } from "./modules/local/bam"

// ///////////////////



// ///////////////
// // gene tagging

// include { htseq } from "./modules/local/tagging"

// include { bam_metrics as count_gene_tags } from "./modules/local/bam"
// count_gene_tags_script = Channel.fromPath("bin/bam/count_gene_tags.py")

// include { bam_metrics as count_reads_per_umi } from "./modules/local/bam"
// count_reads_per_umi_script = Channel.fromPath("bin/bam/reads_per_umi.py")

// include { bam_metrics as count_reads_per_umi_gene } from "./modules/local/bam"
// count_reads_per_umi_gene_script = Channel.fromPath("bin/bam/reads_per_umi_gene.py")

// include { plot_1_arg as plot_gene_tags } from "./modules/local/plot"
// plot_gene_tags_script = Channel.fromPath("bin/plot/gene_tags.py")

// include { bam_filter as bam_filter_gene_tags } from "./modules/local/bam"

// ///////////////

// /////////////////////
// // umis multi mapping

// include { select } from "./modules/local/tagging"
// select_script = Channel.fromPath("bin/select")

// include { duplicates } from "./modules/local/quality_control"
// duplicates_script = Channel.fromPath("bin/duplicates.py")

// include { bam_metrics as count_select } from "./modules/local/bam"
// count_select_script = Channel.fromPath("bin/bam/count_select.py")

// include { plot_1_arg as plot_select } from "./modules/local/plot"
// plot_select_script = Channel.fromPath("bin/plot/select.py")

// include { bam_filter as bam_filter_multimapped_umis } from "./modules/local/bam"
// /////////////////////

// ////////////
// // sequences

// include { bam_metrics as reads_per_barcode_umi } from "./modules/local/bam"
// reads_per_barcode_umi_script = Channel.fromPath("bin/bam/reads_per_barcode_umi.py")

// include { plot_1_arg as plot_balance_barcode } from "./modules/local/plot"
// plot_balance_barcode_script = Channel.fromPath("bin/plot/balance_barcode.py")

// include { plot_1_arg as plot_balance_umi } from "./modules/local/plot"
// plot_balance_umi_script  = Channel.fromPath("bin/plot/balance_umi.py")

// include { plot_1_arg as plot_reads_fraction } from "./modules/local/plot"
// plot_reads_fraction_script = Channel.fromPath("bin/plot/reads_fraction.py")
// ////////////

// //////
// // dge

// include { dge } from "./modules/local/export"
// count_script = Channel.fromPath("bin/count")

// plot_histo_genes_script = Channel.fromPath("bin/plot/histo_genes.py")
// plot_histo_umis_script = Channel.fromPath("bin/plot/histo_umis.py")
// plot_spatial_umis_script  = Channel.fromPath("bin/plot/spatial_umi.py")
// plot_umis_per_barcode_script = Channel.fromPath("bin/plot/umis_per_barcode.py")

// include { plot_1_arg as plot_histo_genes } from "./modules/local/plot"
// include { plot_1_arg as plot_histo_umis } from "./modules/local/plot"
// include { plot_1_arg as plot_umis_per_barcode } from "./modules/local/plot"
// include { plot_2_args as plot_spatial_umis } from "./modules/local/plot"
// //////

// /////////
// // output 

// include { merge_plots } from "./modules/local/export"
// include { rename_coords } from "./modules/local/export"
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
	// PREPARE GENOMES

    if (params.star_index.endsWith(".tar.gz")) {
        ch_star_index = UNTAR_STAR_INDEX ( [ [], params.star_index ] ).untar.map{ row -> [ row[1] ] }
    } else {
        ch_star_index = file(params.star_index)
    }

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

	///////////////////////////////////////////////////////////////////////////
	// ALIGNMENT

	STAR( EXTRACT_BARCODE.out.fastq, ch_star_index )

	///////////////////////////////////////////////////////////////////////////
	// DUPLICATES

	MARK_DUPLICATES( STAR.out.bam )

	///////////////////////////////////////////////////////////////////////////
	// ADD SLIDE-SEQ, ALIGNMENT AND DUPLICATES TAGS

	// slide-seq tags, alignment tag and duplicate tag
	TAG_BAM( MARK_DUPLICATES.out.bam )

	// 3 columns: Matched, Mapped, Reads
	READS_UP_MATCHING( TAG_BAM.out )

	PLOT_UP_ALIGN( READS_UP_MATCHING.out.map{ [ it[0] , it[1], 'up_align' ] } )

	BAM_FILTER_UP_MATCHED( TAG_BAM.out )

	///////////////////////////////////////////////////////////////////////////
	// UMIS PER BARCODE THRESHOLD

	UMIS_PER_BARCODE( BAM_FILTER_UP_MATCHED.out )

	READS_UMIS_PER_BARCODE( UMIS_PER_BARCODE.out )

	READS_UMI_THRESHOLD( UMIS_PER_BARCODE.out )

	PLOT_UMI_THRESHOLD( READS_UMI_THRESHOLD.out.map{ [ it[0] , it[1], params.umi_threshold ] } )

	BAM_FILTER_UMI_THRESHOLD( UMIS_PER_BARCODE.out )

	// ///////////////////////////////////////////////////////////////////////////
	// // HAMMING DISTANCE

	GET_BARCODES( READS_UMIS_PER_BARCODE.out )

	GET_BARCODES.out
		.combine( ch_puck_barcodes )
		.filter{ it[0]["puck"] == it[2] }
		.map{ [ addValue(it[0], "barcodes", it[3]) , it[1] , it[4] ] }
		.set{ ch_hamming_input }
	//ch_hamming_input | view

	HAMMING( ch_hamming_input )

	ch_plot_hamming = HAMMING.out.map{ [ it[0]["name"] , *it ] }.groupTuple()
		.map{
			it[1][0]["barcodes"] == "ordered" ? [it[1][0], *it[2]] : [it[1][1], *it[2].reverse()]
		}
	//ch_plot_hamming | view 

	PLOT_HAMMING_HISTO ( ch_plot_hamming )

	// ///////////////////////////////////////////////////////////////////////////
	// // BARCODE MATCHING

	HAMMING
		.out
		.combine(READS_UMIS_PER_BARCODE.out)
		.filter{ it[0]["name"] == it[2]["name"] }
		.map{ [ it[0] , it[1] , it[3] ] }
		.combine( ch_pucks )
		.filter{ it[0]["puck"] == it[3] }
		.map{ [ *it[0..2] , it[4] ] }
		.set{ ch_matching }
	
	MATCHER( ch_matching )

	PLOT_BARCODE_MATCHING(
		MATCHER.out.metrics
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("barcode_matching") )
	)

	PLOT_HISTO_ERROR (
		MATCHER.out.mapping
			.filter{ it[0]["barcodes"] == "ordered" }
			.combine( Channel.from("histo_errors") )
	)

	MATCHER.out.mapping
	.filter{ it[0]["barcodes"] == "ordered" }
	.combine(BAM_FILTER_UMI_THRESHOLD.out)
	.filter{ it[0]["name"] == it[2]["name"] }
	.map{ [ *it[0..1] , it[3] ] }
	.set{ ch_add_match }
	//ch_add_match | view

	ADD_MATCH( ch_add_match )

	READS_BARCODE_MATCHING( ADD_MATCH.out )

	PLOT_BARCODE_ALIGN(
		reads_barcode_matching
			.out
			.combine( Channel.from("barcode_align") )
	)

	BAM_FILTER_BARCODE_MATCHED( ADD_MATCH.out )

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

