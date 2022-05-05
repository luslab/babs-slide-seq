// ============================================================================
// Program: antisense
// Author: Nourdine Bah <nourdine.bah@crick.ac.uk>
// ============================================================================

#include <tuple>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/gff_io.h>
#include <seqan/bam_io.h>

#include "utils.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct AntisenseOptions
{
	CharString gff;
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
parseCommandLine(AntisenseOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("antisense");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program adds the alignment orientation as a tag");
	setVersion(parser, "1.0");
	setDate(parser, "Apr 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIGFF\\fP \\fIBAM\\fP ");
	addDescription(parser,
		"This programs add the alignment orientation as a tag."
		"The \\fIXF\\fP tag must be set to a gene and nothing else."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "GFF"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_IN"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM_OUT"));
	setValidValues(parser, 0, "gtf GTF gff GFF");
	setValidValues(parser, 1, "bam BAM sam SAM");
	setValidValues(parser, 2, "bam BAM sam SAM");

	// Parse the arguments
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	// Only extract options. The program continues after parseCommandLine()
	if (res != ArgumentParser::PARSE_OK)
	{
		return res;
	}

	// Export arguments and option
	getArgumentValue(options.gff, parser, 0);
	getArgumentValue(options.bam_in, parser, 1);
	getArgumentValue(options.bam_out, parser, 2);

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function GetStrandsFromGFF()
// ----------------------------------------------------------------------------

std::map<std::string, char> GetStrandsFromGFF(CharString gff_path)
{
	GffFileIn gff(toCString(gff_path));

	std::set<seqan::CharString> TYPES = {"gene", "exon", "start_codon", "stop_codon", "five_prime_utr", "three_prime_utr"};

	std::map<std::string, char> strands;
	GffRecord record;
	std::set<seqan::CharString>::iterator it;
	int counter = 1;

	while ( ! atEnd(gff) )
	{
		readRecord(record, gff);

		it = TYPES.find(record.type);

		if ( it == TYPES.end() )
		{
			continue;
		}

		int gene_id_position = GetPosition(record.tagNames, "gene_id");

		if ( -1 != gene_id_position )
		{
			std::string gene_id(toCString(record.tagValues[gene_id_position]));
			strands[gene_id] = record.strand;
		}

		// counter
		if ( counter % (unsigned long)100000 == 0 ) {
			std::cerr << counter << " gff entries" << std::endl;
		}
		counter++;
	}

	close(gff);

	return strands;
}

// ----------------------------------------------------------------------------
// Function AddOrientation()
// ----------------------------------------------------------------------------

void
AddOrientation(CharString bam_path_in, CharString bam_path_out,
	std::map<std::string, char> strands, bool test=false)
{
	// BAM files
	BamFileIn bam_in( toCString(bam_path_in) );
	BamFileOut bam_out(context(bam_in), toCString(bam_path_out));
	BamHeader header;
	readHeader(header, bam_in);
	writeHeader(bam_out, header);
	BamAlignmentRecord record;

	// BAM tags
	CharString gene;
	CharString orientation;

	unsigned long long counter = 0;

	while ( ! atEnd(bam_in) )
	{
		readRecord(record, bam_in);

		BamTagsDict dict(record.tags);

		ExtractTag("XF", dict, gene);
		std::string g(toCString(gene));
		char strand = strands[g];

		if ( '+' == strand && hasFlagRC(record) )
		{
			orientation = "SENSE";
		}

		else if ( '-' == strand && !hasFlagRC(record) )
		{
			orientation = "SENSE";
		}

		else
		{
			orientation = "ANTISENSE";
		}

		appendTagValue(dict, "ao", orientation);
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
	AntisenseOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	// Extract gene strand for GFF
	std::map<std::string, char> strands = GetStrandsFromGFF(options.gff);

	// Add the alignment orientation as a BAM tag
	AddOrientation(options.bam_in, options.bam_out, strands, test);

	std::cerr << "Done." << std::endl;
	return 0;
}


