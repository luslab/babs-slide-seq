// ============================================================================
// Program: count
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
#include "gene.hpp"
#include "counts.hpp"
#include "counter.hpp"

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// This struct stores the options from the command line

struct CountOptions
{
	CharString gff;
	CharString bam;
	CharString dir;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(CountOptions& options, int argc, char const** argv)
{
	// Initiatlize ArgumentParser
	seqan::ArgumentParser parser("count");
	setCategory(parser, "Slide-Seq");
	setShortDescription(parser, "This program produces a DGE matrix from a BAM file");
	setVersion(parser, "1.0");
	setDate(parser, "Apr 2022");

	// Add use and description lines
	addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIGFF\\fP \\fIBAM\\fP ");
	addDescription(parser,
		"Each barcode-UMI pair maps to a gene. "
		"This programs counts the number of times a barcode maps a gene and outputs "
		"a Seurat-style DGE matrix (\\fBfeatures.tsv\\fP, \\fBbarcodes.tsv\\fP, \\fBmatrix.mtv\\fP). "
		"Gene IDs and names are gathered from the input GFF file and the resulting DGE matrix will "
		"contain entries even for undected genes. "
		"The input BAM file must contain the following tags for each record: "
		"\\fIbc\\fP for the barcode, \\fImi\\fP for the UMI and \\fIXF\\fP for the "
		"gene ID. "
		"The \\fIXF\\fP tag must be set to a gene and nothing else. "
		"Each barcode-UMI pair should map to only one gene, so the total counts of the DGE "
		"will be equal to the number of records in the BAM file."
	);

	// Add positional arguments and set their valid file types
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "GFF"));
	addArgument(parser, ArgParseArgument(ArgParseArgument::INPUT_FILE, "BAM"));
	setValidValues(parser, 0, "gtf GTF gff GFF");
	setValidValues(parser, 1, "bam BAM sam SAM");

	// Add the DGE directory name
	addOption(parser, ArgParseOption("d", "directory", "The directory containing the DGE files", ArgParseArgument::STRING));
	setDefaultValue(parser, "d", ".");

	// Parse the arguments
	ArgumentParser::ParseResult res = parse(parser, argc, argv);

	// Only extract options. The program continues after parseCommandLine()
	if (res != ArgumentParser::PARSE_OK)
	{
		return res;
	}

	// Export arguments and option
	getArgumentValue(options.gff, parser, 0);
	getArgumentValue(options.bam, parser, 1);
	getOptionValue(options.dir, parser, "directory");

	return seqan::ArgumentParser::PARSE_OK;
}

// ----------------------------------------------------------------------------
// Function GetGenesFromGFF()
// ----------------------------------------------------------------------------

std::set<Gene> GetGenesFromGFF(CharString gff_path)
{
	GffFileIn gff(toCString(gff_path));

	std::set<seqan::CharString> TYPES = {"gene", "exon", "start_codon", "stop_codon", "five_prime_utr", "three_prime_utr"};

	std::set<Gene> genes;
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
		int gene_name_position = GetPosition(record.tagNames, "gene_name");

		if ( -1 != gene_name_position && -1 != gene_id_position )
		{

			Gene gene = Gene(toCString(record.tagValues[gene_id_position]), toCString(record.tagValues[gene_name_position]));
			genes.insert(gene);
		}

		// counter
		if ( counter % (unsigned long)100000 == 0 ) {
			std::cerr << counter << " gff entries" << std::endl;
		}
		counter++;
	}

	close(gff);

	return genes;
}

// ----------------------------------------------------------------------------
// Function CountFromBAM()
// ----------------------------------------------------------------------------

void CountFromBAM(CharString bam_path, Counter& counter, bool test=false)
{
	BamFileIn bam(toCString(bam_path));
	BamHeader header;
	readHeader(header, bam);
	BamAlignmentRecord rec;

	CharString barcode;
	CharString gene;

	unsigned long long position = 0;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict dict(rec.tags);

		ExtractTag("bc", dict, barcode);
		ExtractTag("XF", dict, gene);

		counter.Increment(toCString(barcode), toCString(gene)); 

		position++;

		if ( (position+1) % (unsigned long long)1000000 == 0 ) {
			std::cerr << "Counting " << (position+1) << "..." << std::endl;
		}

		if (test) { if ( position > 999999 ) { break; } }
	}

	close(bam);
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

// Program entry point

int main(int argc, char const** argv)
{
	// Parse command line
	seqan::ArgumentParser parser;
	CountOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	// If there was an error then the programs exits. The return code is 1 if
	// there were errors and 0 if there were none
	if (res != seqan::ArgumentParser::PARSE_OK)
	{
		return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	bool test = false;

	// Extract gene IDs and symbols for GFF
	std::set<Gene> genes = GetGenesFromGFF(options.gff);

	// Counts for each barcode
	Counter counter = Counter(genes);
	CountFromBAM(options.bam, counter, test);

	// Export DGE matrix
	std::cerr << "Exporting..." << std::endl;
	std::string out_dir(toCString(options.dir));
	std::ofstream feat(out_dir + "/features.tsv");
	std::ofstream bcd(out_dir + "/barcodes.tsv");
	std::ofstream mtx(out_dir + "/matrix.mtx");
	counter.GetDGE(feat, bcd, mtx);
	feat.close();
	bcd.close();
	mtx.close();

	std::cerr << "Done." << std::endl;
	return 0;
}

