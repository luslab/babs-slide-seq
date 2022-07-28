// ============================================================================
// Program: umis_per_barcode
// Author: Nourdine Bah <nourdine.bah@crick.ac.uk>
// ============================================================================

#include <tuple>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

#include "utils.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct UmisPerBarcodeOptions
{
	CharString bam_in;
	CharString bam_out;
	int threshold;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(UmisPerBarcodeOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("umis_per_barcode");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program counts the number of UMIs per barcode");
	setVersion(parser, "1.0");
	setDate(parser, "Mar 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM_IN\\fP \\fIBAM_OUT\\fP ");
	addDescription(parser,
		"A bead barcode has multiple UMIs. "
		"The number of UMIs per barcode is a useful metrics. "
		"This programs counts the number of UMIs per barcode and adds this number "
		"as a BAM tag named \\fIbn\\fP. "
		"The programs also adds a \\fIbt\\fP BAM tag that can either be \\fBPASS\\fP "
		"or \\fBTOO_LOW\\fP based on the \\fBthreshold\\fP option."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_IN"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "BAM_OUT"));
	setValidValues(parser, 0, "bam BAM sam SAM");
	setValidValues(parser, 1, "bam BAM sam SAM");

	// Add the threshold option
	addOption(parser, ArgParseOption("t", "threshold", "The mininum number of UMIs per barcode", ArgParseArgument::INTEGER));
	setDefaultValue(parser, "t", 0);

	// Parse the arguments
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	// Only extract options. The program continues after parseCommandLine()
	if (res != ArgumentParser::PARSE_OK)
	{
		return res;
	}

	// Export arguments
	getArgumentValue(options.bam_in, parser, 0);
	getArgumentValue(options.bam_out, parser, 1);
	getOptionValue(options.threshold, parser, "threshold");

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function GetUMIsFromBam()
// ----------------------------------------------------------------------------

std::map<CharString, std::set<CharString>>
GetUMIsFromBAM(CharString bam_path, bool test=false)
{
	// BAM tags values
	CharString up_status;
	CharString align_status;
	CharString barcode;
	CharString umi;
	CharString umi_status;

	// A set of UMIs for each barcode
	std::map<CharString, std::set<CharString>> barcodes;
	//std::map<CharString, CharString> status;

	// BAM file
	BamFileIn bam( toCString(bam_path) );
	BamHeader header;
	BamAlignmentRecord rec;
	readHeader(header, bam);

	unsigned long long counter = 0;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict dict(rec.tags);

		ExtractTag("us", dict, up_status);
		ExtractTag("as", dict, align_status);

		if ( "MATCHED" == up_status && "MAPPED" == align_status )
		{
			ExtractTag("bc", dict, barcode);
			ExtractTag("mi", dict, umi);

			barcodes[barcode].insert(umi);
		}

		counter++;

		if ( (counter+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Extracting " << (counter+1) << "..." << std::endl;
		}

		if (test) { if ( counter > 999999 ) { break; } }
	}

	close(bam);

	return barcodes;
}

// ----------------------------------------------------------------------------
// Function AddTags()
// ----------------------------------------------------------------------------

void
AddTags(CharString bam_in_path, CharString bam_out_path,
	std::map<CharString, CharString>& status, std::map<CharString, unsigned long long>& counts,
	bool test=false)
{
	// BAM tags values
	CharString up_status;
	CharString align_status;
	CharString barcode;

	// BAM files
	BamFileIn bam_in( toCString(bam_in_path) );
	BamFileOut bam_out(context(bam_in), toCString(bam_out_path));

	// Header
	BamHeader header;
	readHeader(header, bam_in);
	writeHeader(bam_out, header);

	BamAlignmentRecord rec;
	unsigned long long counter = 0;

	while ( ! atEnd(bam_in) )
	{
		readRecord(rec, bam_in);
		BamTagsDict dict(rec.tags);

		ExtractTag("us", dict, up_status);
		ExtractTag("as", dict, align_status);

		if ( "MATCHED" == up_status && "MAPPED" == align_status )
		{
			ExtractTag("bc", dict, barcode);
			appendTagValue(dict, "bn", std::to_string(counts[barcode]));
			appendTagValue(dict, "bt", status[barcode]);
		}

		else
		{
			appendTagValue(dict, "bn", "NULL");
			appendTagValue(dict, "bt", "NULL");
		}

		//tagsToBamRecord(rec, dict);
		writeRecord(bam_out, rec);

		counter++;

		if ( (counter+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Tagging " << (counter+1) << "..." << std::endl;
		}

		if (test) { if ( counter > 999999 ) { break; } }
	}

	close(bam_in);
	close(bam_out);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

// Program entry point

int main(int argc, char const** argv)
{
	// Parse command line
	seqan::ArgumentParser parser;
	UmisPerBarcodeOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	// Get all the UMIs for each barcode
	std::map<CharString, std::set<CharString>> barcodes = GetUMIsFromBAM(options.bam_in, test);

	// Counts the number of UMIs and set the status tag
	// (PASS or TWO_LOW based on the threshold)
	std::map<CharString, CharString> status;
	std::map<CharString, unsigned long long> counts;
	for (auto& [barcode, umis] : barcodes)
	{
		status[barcode] = (umis.size() < options.threshold) ? "TOO_LOW" : "PASS";
		counts[barcode] = umis.size();
	}

	// Add tags to BAM files
	AddTags(options.bam_in, options.bam_out, status, counts, test);

	return 0;
}

