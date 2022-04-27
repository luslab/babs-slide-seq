// ============================================================================
// Program: tag_bam
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

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct TagBamOptions
{
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
parseCommandLine(TagBamOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("tag_bam");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program converts record names to bam tags");
	setVersion(parser, "1.0");
	setDate(parser, "Apr 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM_IN\\fP \\fIBAM_OUT\\fP ");
	addDescription(parser,
		"The information (barcode, UMI...) is contained in the record name. "
		"The program parses the record name and converts the information to BAM tags."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_IN"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUT_FILE, "BAM_OUT"));
	setValidValues(parser, 0, "bam BAM sam SAM");
	setValidValues(parser, 1, "bam BAM sam SAM");

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

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function AddTags()
// ----------------------------------------------------------------------------

void AddTags(CharString bam_in_path, CharString bam_out_path, bool test=false)
{
	// BAM files
	BamFileIn bam_in(toCString(bam_in_path));
	BamFileOut bam_out(context(bam_in), toCString(bam_out_path));
	BamHeader header;
	readHeader(header, bam_in);
	writeHeader(bam_out, header);
	BamAlignmentRecord record;

	// The parsing regex
	std::string regex;
	regex.append("([^:]+)"); // length_status
	regex.append(":");
	regex.append("([A-Z]+)"); // up
	regex.append(":");
	regex.append("(-*[0-9]+)"); // up_hamming
	regex.append(":");
	regex.append("([^:]+)"); // up_status
	regex.append(":");
	regex.append("([A-Z]+)"); // barcode
	regex.append(":");
	regex.append("([A-Z]+)"); // umi
	regex.append(":");
	regex.append("(.*)$"); // read_name
	std::regex rgx(regex);

	// The information
	std::string length_status;
	std::string up;
	int up_hamming;
	std::string up_status;
	std::string barcode;
	std::string umi;
	std::string read_name;
	std::string mapping_status;
	std::string duplicate_status;

	unsigned long long counter = 0;
	std::smatch m;

	while ( ! atEnd(bam_in) )
	{
		// Current record
		readRecord(record, bam_in);

		// Parse the name and gets the information
		std::string qname( toCString(record.qName) );
		if ( std::regex_match(qname, m, rgx) )
		{
			length_status = m[1];
			up = m[2];
			up_hamming = std::stoi(m[3]);
			up_status = m[4];
			barcode = m[5];
			umi = m[6];
			read_name = m[7];
		}
		else
		{
			std::cerr << "Unable to match regex" << std::endl;
			std::cerr << regex << std::endl;
			std::cerr << qname << std::endl;
			exit(EXIT_FAILURE);
		}

		// Remove the information from the record name
		record.qName = read_name.c_str();

		// Mapping and duplicates tags
		if ( hasFlagUnmapped(record) )
		{
			mapping_status = "UNMAPPED";
			duplicate_status = "NULL";
		}
		else
		{
			mapping_status = "MAPPED";

			if ( hasFlagDuplicate(record) )
			{
				duplicate_status = "DUPLICATE";
			}
			else
			{
				duplicate_status = "PRIMARY";
			}
		}

		// Add the tags to the BAM record
		BamTagsDict tagsDict(record.tags);
		appendTagValue(tagsDict, "ls", length_status.c_str());
		appendTagValue(tagsDict, "up", up.c_str());
		appendTagValue(tagsDict, "uh", up_hamming);
		appendTagValue(tagsDict, "us", up_status.c_str());
		appendTagValue(tagsDict, "bc", barcode.c_str());
		appendTagValue(tagsDict, "mi", umi.c_str());
		appendTagValue(tagsDict, "as", mapping_status.c_str());
		appendTagValue(tagsDict, "dp", duplicate_status.c_str());
		tagsToBamRecord(record, tagsDict);

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
	TagBamOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	AddTags(options.bam_in, options.bam_out, test);

	return 0;
}

