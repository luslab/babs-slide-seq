// ============================================================================
// Program: add_match
// Author: Nourdine Bah <nourdine.bah@crick.ac.uk>
// ============================================================================

#include <tuple>
#include <map>
#include <set>
#include <fstream>
#include <vector>
#include <regex>
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

struct AddMatchOptions
{
	CharString csv;
	CharString bam_in;
	CharString bam_out;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AddMatchOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("add_match");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program adds the matched puck barcode as a BAM tag");
	setVersion(parser, "1.0");
	setDate(parser, "Apr 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fICSV\\fP \\fIBAM_IN\\fP \\fIBAM_OUT\\fP ");
	addDescription(parser,
		"Sequencing barcodes are matched with puck barcodes. "
		"This program takes the result of the matching and add the puck barcode to the corresponding "
		"sequencing barcode. "
		"The puck barcode sequence is added with the \\fIbm\\fP tag. "
		"A BAM tag named \\fIbs\\fP is also added. "
		"Its value is says if the sequencing barcode could be matched and is either "
		"\\fBMATCHED\\fP or \\fBUNMATCHED\\fP. "
		"The CSV file must contain two columns with a header. "
		"The first column is the sequencing barcode and the second is the puck barcode."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "CSV"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_IN"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "BAM_OUT"));
	setValidValues(parser, 0, "csv CSV");
	setValidValues(parser, 1, "bam BAM sam SAM");
	setValidValues(parser, 2, "bam BAM sam SAM");

	// Parse the arguments
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	// Only extract options. The program continues after parseCommandLine()
	if (res != ArgumentParser::PARSE_OK)
	{
		return res;
	}

	// Export arguments
	getArgumentValue(options.csv, parser, 0);
	getArgumentValue(options.bam_in, parser, 1);
	getArgumentValue(options.bam_out, parser, 2);

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function GetMatchingsFromCSV()
// ----------------------------------------------------------------------------

std::tuple<CharString, std::map<CharString, CharString>>
GetMatchingsFromCSV(CharString csv_path)
{
	std::map<CharString, CharString> mapping;

	std::ifstream csv( toCString(csv_path) );
	std::string line;
	CharString seq_barcode, puck_barcode;

	std::set<int> lengths;

	int i = 0;

	// The parsing regex
	std::regex regex("^([A-Z]+),([A-Z]+).*$");
	std::smatch m;

	while( std::getline(csv, line) )
	{
		// header
		if ( i == 0 )
		{
			i++;
			continue;
		}

		if ( std::regex_match(line, m, regex) )
		{
			std::string s_bcd(m[1]);
			std::string p_bcd(m[2]);
			seq_barcode = s_bcd;
			puck_barcode = p_bcd;
		}
		else
		{
			std::cerr << "Error: cannot parse line " << i << std::endl;
			exit(EXIT_FAILURE);
		}

		if ( length(seq_barcode) != length(puck_barcode) )
		{
			std::cerr << "Error: the barcodes pair " << i << " have different lengths" << std::endl;
			exit(EXIT_FAILURE);
		}
		else
		{
			lengths.insert( length(seq_barcode) );
			mapping[seq_barcode] = puck_barcode;
		}
	}

	csv.close();


	int length;
	std::set<int>::iterator lit;

	if ( 1 != lengths.size() )
	{
		std::cerr << "Error: not all the barcodes pair have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	else
	{
		lit = lengths.begin();
		length = *lit;
	}

	return std::make_tuple(DummySeq(length), mapping);
}

// ----------------------------------------------------------------------------
// Function AddTag()
// ----------------------------------------------------------------------------

void
AddTag(CharString bam_path_in, CharString bam_path_out, CharString dummy,
	std::map<CharString, CharString> mapping, bool test=false)
{
	// BAM files
	BamFileIn bam_in( toCString(bam_path_in) );
	BamFileOut bam_out(context(bam_in), toCString(bam_path_out));
	BamHeader header;
	readHeader(header, bam_in);
	writeHeader(bam_out, header);
	BamAlignmentRecord record;

	// BAM tags
	CharString barcode;
	CharString up_status;

	CharString puck_seq;
	std::string status;

	unsigned long long counter = 0;

	while ( ! atEnd(bam_in) )
	{
		readRecord(record, bam_in);

		BamTagsDict dict(record.tags);

		ExtractTag("bc", dict, barcode);
		ExtractTag("us", dict, up_status);

		if ( up_status == "MATCHED" )
		{
			if ( mapping.find(barcode) != mapping.end() )
			{
				status = "MATCHED";
				puck_seq = mapping[barcode];
			}
			else
			{
				status = "UNMATCHED";
				puck_seq = dummy;
			}
		}
		else
		{
				status = "NULL";
				puck_seq = dummy;
		}

		appendTagValue(dict, "bs", status.c_str());
		appendTagValue(dict, "bm", toCString(puck_seq));
		tagsToBamRecord(record, dict);
		writeRecord(bam_out, record);

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
	AddMatchOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	// Extract the matching from the CSV file
	std::tuple<CharString, std::map<CharString, CharString>> matchings = GetMatchingsFromCSV(options.csv);
	CharString dummy = std::get<0>(matchings);
	std::map<CharString, CharString> mapping = std::get<1>(matchings);

	// Add the puck sequence as a BAM tag
	AddTag(options.bam_in, options.bam_out, dummy, mapping, test);

	return 0;
}

