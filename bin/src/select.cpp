// ============================================================================
// Program: select
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
#include "record.hpp"
#include "molecule.hpp"
#include "molecules.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct SelectOptions
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
parseCommandLine(SelectOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("select");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program tries to match a barcode-UMI pair to a unique gene");
	setVersion(parser, "1.0");
	setDate(parser, "Mar 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIBAM_IN\\fP \\fIBAM_OUT\\fP ");
	addDescription(parser,
		"A barcode-UMI pair is usually associated with multiple alignment records, "
		"and each record is associated with a gene. "
		"So, a barcode-UMI pair can be associated with multiple genes. "
		"This programs tries to match a unique gene to a barcode-UMI pair. "
		"If the pair is associated with a unique gene, then it takes that gene. "
		"If the pair has got several genes, then it searches for a maximum alignment "
		"score among the records and if it can find a maximum then it takes the "
		"record associated with it and discards the other ones. "
		"If it can't find a maximum, then it discards all the records. "
		"The program takes a BAM file a an input and outputs the selected records "
		"in another BAM file. "
		"\\fBThe input BAM file must not contain alignment records that don't map "
		"any gene\\fP. "
		"The input BAM file must contain the following tags for each record: "
		"\\fIbc\\fP for the barcode, \\fImi\\fP for the UMI, \\fIXF\\fP for the "
		"gene and \\fIAS\\fP for the alignment score. "
		"The \\fIXF\\fP tag must be set to a gene and nothing else."
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
// Function ExtractMoleculesFromBAM()
// ----------------------------------------------------------------------------

Molecules ExtractMoleculesFromBAM(CharString bam_path, bool test=false)
{
	unsigned long long position = 0;

	BamFileIn bam(toCString(bam_path));
	BamHeader header;
	readHeader(header, bam);

	BamAlignmentRecord rec;

	Molecules molecules;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict dict(rec.tags);

		CharString barcode;
		CharString umi;
		CharString gene;
		long score;

		ExtractTag("bc", dict, barcode);
		ExtractTag("mi", dict, umi);
		ExtractTag("XF", dict, gene);
		ExtractTag("AS", dict, score);

		Record record = Record(position, toCString(rec.qName), score, toCString(gene));
		Molecule molecule = Molecule(toCString(barcode), toCString(umi), record);
		molecules.Insert(molecule);

		position++;

		if ( (position+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Extracting " << (position+1) << "..." << std::endl;
		}

		if (test) { if ( position > 9999999 ) { break; } }
	}

	close(bam);

	return molecules;
}

// ----------------------------------------------------------------------------
// Function WriteTaggedRecords()
// ----------------------------------------------------------------------------

void WriteTaggedRecords(
		CharString bam_path_in, CharString bam_path_out,
		std::map<unsigned long long, std::string>& values, char* tag, bool test=false)
{
	BamFileIn bam_in(toCString(bam_path_in));
	BamFileOut bam_out(context(bam_in), toCString(bam_path_out));

	BamHeader header;
	readHeader(header, bam_in);
	writeHeader(bam_out, header);

	BamAlignmentRecord rec;

	unsigned long long position = 0;

	while ( ! atEnd(bam_in) )
	{
		readRecord(rec, bam_in);
		BamTagsDict dict(rec.tags);
		appendTagValue(dict, tag, values[position]);
		tagsToBamRecord(rec, dict);
		writeRecord(bam_out, rec);

		position++;

		if ( (position+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Writing " << (position+1) << "..." << std::endl;
		}

		if (test) { if ( position > 9999999 ) { break; } }
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
	SelectOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	// The Molecule class allows to try to find a gene for a barcode-UMI pair
	Molecules molecules = ExtractMoleculesFromBAM(options.bam_in, test);

	// For each alignment record, assign a mapping tag:
	// UNIQUE, INCLUDED, EXCLUDED, UNRESOLVED
	std::map<unsigned long long, std::string> tags;
	unsigned long long i = 0;
	for (auto& molecule : molecules)
	{
		Molecule mol = Molecule(molecule);
		mol.ComputeFrequencies();

		for (auto& [pos, tag] : mol.GetFrequencyBasedRecordTags())
		{
			tags[pos] = tag;
		}

		i++;
		if ( (i+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Processing " << (i+1) << "..." << std::endl;
		}
	}

	// Write the records with the mapping tags in a new BAM file
	WriteTaggedRecords(options.bam_in, options.bam_out, tags, "cs", test);

	return 0;
}

