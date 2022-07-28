// ============================================================================
// Program: extract_barcode
// Author: Nourdine Bah <nourdine.bah@crick.ac.uk>
// ============================================================================

#include <map>
#include <set>
#include <fstream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "utils.hpp"
#include "read_structure.hpp"
#include "read1.hpp"
#include "read2.hpp"
#include "transcript.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct ExtractBarcodeOptions
{
	CharString fastq1;
	CharString fastq2;
	CharString fastq;
	CharString sample;
	CharString structure;
	int up_distance;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(ExtractBarcodeOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("extract_barcode");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program extracts bead barcodes, UMIs and transcripts");
	setVersion(parser, "1.0");
	setDate(parser, "Mar 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIFASTQ1\\fP \\fIFASTQ2\\fP ");
	addDescription(parser,
		"This programs takes Read 1 and Read 2 as input. It uses the read structure "
		"definition to extract bead barcodes and UMIs from Read 1. Then, it extracts the "
		"transcript from Read 2 and creates a FASTQ whose reads are the transcripts "
		"(Read 2) and descriptions contains bead barcodes, UMIs and other information. "
		"You need to specify a sample name, a read structure definition and a maximum "
		"hamming distance for the UP primer."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READ1"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "READ2"));
	setValidValues(parser, 0, "fq fq.gz fastq fastq.gz");
	setValidValues(parser, 1, "fq fq.gz fastq fastq.gz");

	// Add section with sample information
	addSection(parser, "Sample");
	addOption(parser, ArgParseOption("s", "sample", "The sample name", ArgParseArgument::STRING));
	setDefaultValue(parser, "s", "Sample");
	addOption(parser, ArgParseOption("f", "fastq", "The path of the output FASTQ", ArgParseArgument::OUTPUT_FILE));
	setDefaultValue(parser, "f", "Sample.fastq.gz");

	// Add section with parameters
	addSection(parser, "Parameters");
	addOption(parser, ArgParseOption("r", "read-structure", "Read structure definition", ArgParseArgument::STRING));
	setDefaultValue(parser, "r", "8C18U6C3X8M");
	addOption(parser, ArgParseOption("d", "max-distance", "Maximum hamming distance for UP primer", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "d", 3);

	// Parse the arguments
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	// Only extract options. The program continues after parseCommandLine()
	if (res != ArgumentParser::PARSE_OK)
	{
		return res;
	}

	// Export arguments
	getArgumentValue(options.fastq1, parser, 0);
	getArgumentValue(options.fastq2, parser, 1);

	// Export options
	getOptionValue(options.sample, parser, "sample");
	getOptionValue(options.fastq, parser, "fastq");
	getOptionValue(options.structure, parser, "read-structure");
	getOptionValue(options.up_distance, parser, "max-distance");

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function Extract()
// ----------------------------------------------------------------------------

void Extract(
		SeqFileIn& fastq1, SeqFileIn& fastq2, ReadStructure structure, SeqFileOut& fastq,
		CharString sample, int max_distance, std::map<std::string, unsigned long long>& metrics,
		std::map<int, unsigned long long>& distances)
{
	// Read name, sequence and quality
	CharString meta1;
	CharString meta2;
	CharString sequence1;
	CharString sequence2;
	CharString qual1;
	CharString qual2;

	// Iterate over both FASTQ files
	while ( ! atEnd(fastq1) )
	{
		// Extract read name, sequence and quality
		readRecord(meta1, sequence1, qual1, fastq1);
		readRecord(meta2, sequence2, qual2, fastq2);

		// Parse Read1
		Read1 read1 = Read1(structure, toCString(sequence1), max_distance);

		// The total number of reads in the FASTQ file
		metrics["Total reads"]++;

		// Record an histogram of the UP primer hamming distances from the original
		// sequence
		distances[read1.GetUPDist()]++;

		// Discard the reads that are not long enough and have an UP primer whose
		// the hamming distance from the original sequence exeeds the threshold
		// value
		if ( "LONG_ENOUGH" == read1.GetSeqStatus() )
		{
			metrics["Long enough reads"]++;

			// Create a Transcript object
			Read2 read2 = Read2(toCString(meta2), toCString(sequence2), toCString(qual2));
			Transcript transcript = Transcript(read1, read2);

			if ( "MATCHED" == read1.GetUPStatus() )
			{ 
				metrics["UP primer match"]++;

			}
			else
			{
				metrics["UP primer non match"]++;
			}

			// Export the transcript if the Read 1 is long enough and matches
				// the UP primer sequence
			writeRecord(fastq, transcript.GetMeta(), transcript.GetSeq(), transcript.GetQual());
		}

		else {
			metrics["Too short reads"]++;
		}

		// Display the counter of the number of reads
		if ( metrics["Total reads"] % (unsigned long long)1000000 == 0 )
		{
			std::cerr << metrics["Total reads"] << std::endl;
		}
	}
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

// Program entry point

int main(int argc, char const** argv)
{
	// Parse command line
	seqan::ArgumentParser parser;
	ExtractBarcodeOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	// Create a ReadStructure object from the definition
	ReadStructure structure = ReadStructure(toCString(options.structure));

	// Read1 and Read2 FASTQ files and the output FASTQ file
	SeqFileIn fastq1(toCString(options.fastq1));
	SeqFileIn fastq2(toCString(options.fastq2));
	SeqFileOut fastq(toCString(options.fastq));

	// The metrics and the histogram of the UP primer distances
	std::map<std::string, unsigned long long> metrics;
	std::map<int, unsigned long long> distances;

	// Extract the bead barcodes, UP primer and UMIs
	Extract(fastq1, fastq2, structure, fastq, options.sample, options.up_distance, metrics, distances);

	// Close the FASTQ files
	close(fastq1);
	close(fastq2);
	close(fastq);

	// Sample name to use later
	std::string sample = std::string(toCString(options.sample));

	// Export metrics
	metrics["UP primer length"] = structure.GetLength('U');
	metrics["Bead barcode length"] = structure.GetLength('C');
	metrics["UMI length"] = structure.GetLength('M');
	WriteCounter(sample, sample + ".extract_barcode.csv", metrics);

	// Export histogram of hamming distances
	WriteCounter(sample, sample + ".extract_barcode.up_distances.csv", distances);

	return 0;
}

