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

include { HAMMING            } from "./modules/local/integration"
include { PLOT_HAMMING_HISTO } from "./modules/local/plot"

// ///////////////////

include { MATCHER                                  } from "./modules/local/integration"
include { PLOT as PLOT_BARCODE_MATCHING            } from "./modules/local/plot"
include { PLOT as PLOT_HISTO_ERROR                 } from "./modules/local/plot"
include { ADD_MATCH                                } from "./modules/local/integration"
include { BAM_METRICS as READS_BARCODE_MATCHING    } from "./modules/local/bam"
include { PLOT as PLOT_BARCODE_ALIGN               } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_BARCODE_MATCHED } from "./modules/local/bam"

// ///////////////////
// ///////////////
// // gene tagging

include { HTSEQ                                   } from "./modules/local/tagging"
include { BAM_METRICS as COUNT_GENE_TAGS          } from "./modules/local/bam"
include { PLOT as PLOT_GENE_TAGS                  } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_GENE_TAGS      } from "./modules/local/bam"
include { BAM_METRICS as COUNT_READS_PER_UMI      } from "./modules/local/bam"
include { BAM_METRICS as COUNT_READS_PER_UMI_GENE } from "./modules/local/bam"

// ///////////////
// /////////////////////
// // umis multi mapping

include { SELECT                                    } from "./modules/local/tagging"
include { DUPLICATES                                } from "./modules/local/quality_control"
include { BAM_METRICS as BAM_METRICS_COUNT_SELECT   } from "./modules/local/bam"
include { PLOT as PLOT_SELECT                       } from "./modules/local/plot"
include { BAM_FILTER as BAM_FILTER_MULTIMAPPED_UMIS } from "./modules/local/bam"

// /////////////////////
// ////////////
// // sequences

include { BAM_METRICS as BAM_METRICS_READS_PER_BARCODE_UMI } from "./modules/local/bam"
include { PLOT as PLOT_BALANCE_BARCODE                     } from "./modules/local/plot"
include { PLOT as PLOT_BALANCE_UMI                         } from "./modules/local/plot"
include { PLOT as PLOT_READS_FRACTION                      } from "./modules/local/plot"

// ////////////
// //////
// // dge

include { DGE                           } from "./modules/local/export"
include { PLOT as PLOT_HISTO_GENES      } from "./modules/local/plot"
include { PLOT as PLOT_HISTO_UMIS       } from "./modules/local/plot"
include { PLOT as PLOT_UMIS_PER_BARCODE } from "./modules/local/plot"
include { PLOT as PLOT_SPATIAL_UMIS     } from "./modules/local/plot"

// //////
// /////////
// // output 

include { MERGE_PLOTS   } from "./modules/local/export"
include { RENAME_COORDS } from "./modules/local/export"

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
		READS_BARCODE_MATCHING
			.out
			.combine( Channel.from("barcode_align") )
	)

	BAM_FILTER_BARCODE_MATCHED( ADD_MATCH.out )

	// ///////////////////////////////////////////////////////////////////////////
	// // GENE COUNT

	HTSEQ( BAM_FILTER_BARCODE_MATCHED.out, params.gtf )

	COUNT_GENE_TAGS( HTSEQ.out.bam )

	PLOT_GENE_TAGS( COUNT_GENE_TAGS.out.map{ [ it[0] , it[1], '' ] } )

	BAM_FILTER_GENE_TAGS( HTSEQ.out.bam )

	COUNT_READS_PER_UMI( BAM_FILTER_GENE_TAGS.out )

	COUNT_READS_PER_UMI_GENE( BAM_FILTER_GENE_TAGS.out )

	// ///////////////////////////////////////////////////////////////////////////
	// // UMIS MAPPINGS

	SELECT( BAM_FILTER_GENE_TAGS.out.map{ it[0..1] } )

	SELECT.out
		.unique_reads
		.concat(SELECT.out.resolved_reads)
		.concat(SELECT.out.unresolved_reads)
		.combine(MERGE_LANES.out)
		.filter{ it[0]["name"] == it[3]["name"] }
		.map{ [ addValue(it[0], "status", it[1]) , it[2] , *it[4..5] ] }
		.set{ ch_duplicates }

	DUPLICATES( ch_duplicates )

	BAM_METRICS_COUNT_SELECT( SELECT.out.bam )

	PLOT_SELECT ( BAM_METRICS_COUNT_SELECT.out.map{ [ it[0] , it[1], '' ] } )

	BAM_FILTER_MULTIMAPPED_UMIS( SELECT.out.bam )

	// ///////////////////////////////////////////////////////////////////////////
	// // SEQUENCES

	BAM_METRICS_READS_PER_BARCODE_UMI( BAM_FILTER_MULTIMAPPED_UMIS.out )

	PLOT_BALANCE_BARCODE ( BAM_METRICS_READS_PER_BARCODE_UMI.out.map{ [ it[0] , it[1], '' ] } )

	PLOT_BALANCE_UMI ( BAM_METRICS_READS_PER_BARCODE_UMI.out.map{ [ it[0] , it[1], '' ] } )

	PLOT_READS_FRACTION ( BAM_METRICS_READS_PER_BARCODE_UMI.out.map{ [ it[0] , it[1], '' ] } )

	// ///////////////////////////////////////////////////////////////////////////
	// // EXPRESSION MATRIX

	DGE( BAM_FILTER_MULTIMAPPED_UMIS.out.map{it[0..1]}, params.gtf )

	PLOT_HISTO_UMIS ( DGE.out.map{ [ it[0] , it[1], '' ] } )

	PLOT_HISTO_GENES ( DGE.out.map{ [ it[0] , it[1], '' ] } )

	PLOT_UMIS_PER_BARCODE ( DGE.out.map{ [ it[0] , it[1], '' ] } )

	// ///////////////
	// // spatial umis
	PLOT_SPATIAL_UMIS(
		DGE.out
			.combine(MATCHER.out.coords)
			.filter{ it[2]["barcodes"] == "ordered" }
			.filter{ it[0]["name"] == it[2]["name"] }
			.map{ [ * it[0..1] , it[3] ] }
			.combine( Channel.from("spatial_umis") )
	)
	// ////////////////

	// ///////////////////////////////////////////////////////////////////////////
	// // OUTPUT
	
	PLOT_BARCODE_EXTRACTION
		.out
		.pdf
		.concat(
			PLOT_UP_MATCHING.out.pdf,
			PLOT_UP_ALIGN.out.pdf,
			PLOT_UMI_THRESHOLD.out.pdf,
			PLOT_BARCODE_ALIGN.out.pdf,
			PLOT_GENE_TAGS.out.pdf,
			PLOT_SELECT.out.pdf,
			PLOT_BALANCE_BARCODE.out.pdf,
			PLOT_BALANCE_UMI.out.pdf,
			PLOT_READS_FRACTION.out.pdf,
			PLOT_HISTO_UMIS.out.pdf,
			PLOT_HISTO_GENES.out.pdf,
			PLOT_BARCODE_MATCHING.out.pdf,
			PLOT_HISTO_ERROR.out.pdf,
			PLOT_HAMMING_HISTO.out.pdf,
			PLOT_UMIS_PER_BARCODE.out.pdf,
			PLOT_SPATIAL_UMIS.out.pdf
		)
		.map{ [ it[0]["name"] , it[1] ] }
		.groupTuple()
		.set{ ch_merged_plots }
	
	MERGE_PLOTS(ch_merged_plots)

	RENAME_COORDS( MATCHER.out.coords.filter{ it[0]["barcodes"] == "ordered"} )
}

