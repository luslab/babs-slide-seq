#include <map>
#include <set>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>
#include <seqan/store.h>

using namespace seqan;

typedef StringSet<CharString> MyCargo;
typedef IntervalAndCargo<int, MyCargo> MyInterval;
typedef std::map<CharString, String<MyInterval>> Intervals;
typedef IntervalTree<int, MyCargo> MyTree;
typedef std::map<CharString, MyTree> Trees;

///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int get_position(StringSet<CharString> tagNames, CharString tagName)
{
	int position = -1;

	for (int i=0; i<length(tagNames); i++)
	{
		if ( tagName == tagNames[i] )
		{
			position = i;
		}
	}

	return position;
}

///////////////////////////////////////////////////////////////////////////////
enum FUNCTION
{
	INTERGENIC	= 1 << 1,
	GENIC			= 1 << 2,
	UTR			= 1 << 3,
	CODING		= 1 << 4,
	IGNORED		= 1 << 5,
	UNKNOWN		= 1 << 6
};

// 4	: INTRONIC
// 20	: INTRONIC & CODING
// 28 : INTRONIC & CODING & UTR

/////////////////
int functionFlag(
	std::set<CharString> functions,
	std::map<CharString, int> flags
)
{
	int flag = 0;

	std::set<CharString>::iterator it;
	std::map<CharString, int>::iterator iit;

	for ( it = functions.begin(); it != functions.end(); ++it)
	{
		iit = flags.find(*it);
		if ( iit == flags.end() )
		{
			flag = flag | FUNCTION::UNKNOWN;
		}
		else
		{
			flag = flag | flags[*it];
		}
	}

	return flag;
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	std::map<CharString, int> FLAGS;
	FLAGS["CDS"]					= FUNCTION::IGNORED;
	FLAGS["Selenocysteine"]		= FUNCTION::IGNORED;
	FLAGS["exon"]					= FUNCTION::CODING;
	FLAGS["five_prime_utr"]		= FUNCTION::UTR;
	FLAGS["gene"]					= FUNCTION::GENIC;
	FLAGS["start_codon"]			= FUNCTION::CODING;
	FLAGS["stop_codon"]			= FUNCTION::CODING;
	FLAGS["three_prime_utr"]	= FUNCTION::UTR;
	FLAGS["transcript"]			= FUNCTION::IGNORED;

	////////////////////////////////////////////////////////////////////////////
	GffFileIn gff(argv[1]);
	GffRecord record;

	Intervals map;
	Intervals::iterator i_it;

	////////////////////////////////////////////////////////////////////////////
	std::set<CharString> TYPES = {
		"gene", "exon",
		"start_codon", "stop_codon",
		"five_prime_utr", "three_prime_utr"
	};

	std::set<CharString>::iterator t_it;

	////////////////////////////////////////////////////////////////////////////
	
	int g = 1;

	while ( ! atEnd(gff) )
	{
		readRecord(record, gff);

		/////////////////////////////////////////////////////////////////////////
		t_it = TYPES.find(record.type);
		if ( t_it == TYPES.end() )
		{
			continue;
		}

		/////////////////////////////////////////////////////////////////////////
		i_it = map.find(record.ref);

		int namePos = get_position(record.tagNames, "gene_name");
		int idPos = get_position(record.tagNames, "gene_id");

		if ( -1 != namePos && -1 != idPos )
		{
			MyCargo cargo;
			resize(cargo, 4);
			assignValue(cargo, 0, record.type);
			assignValue(cargo, 1, (CharString)record.strand);
			assignValue(cargo, 2, record.tagValues[namePos]);
			assignValue(cargo, 3, record.tagValues[idPos]);

			if ( i_it == map.end() )
			{
				String<MyInterval> intervals;
				MyInterval interval(record.beginPos, record.endPos, cargo);
				appendValue(intervals, interval);
				map[record.ref] = intervals;

			}
			else
			{
				MyInterval interval(record.beginPos, record.endPos, cargo);
				appendValue(map[record.ref], interval);
			}
		}
		// counter
		if ( g % (unsigned long)100000 == 0 ) {
			std::cerr << g << " gff entries" << std::endl;
		}
		g++;
	}

	close(gff);

	Trees trees;
	for (i_it = map.begin(); i_it != map.end(); ++i_it)
	{
		IntervalTree<int, MyCargo> tree(i_it->second);
		trees[i_it->first] = tree;
	}
	////////////////////////////////////////////////////////////////////////////

	BamFileIn bamFileIn(argv[2]);
	BamFileOut bamFileOut(context(bamFileIn), argv[3]);

	BamHeader header;
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	BamAlignmentRecord rec;

	std::set<CharString> types;
	std::set<CharString> symbols;
	std::set<CharString> ids;

	std::set<CharString>::iterator it;

	int n_genes;
	std::string symbol;
	std::string id;
	std::string function;
	std::string status;

	unsigned long counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(rec, bamFileIn);
		BamTagsDict tagsDict(rec.tags);

		// MAPPED
		if ( ! hasFlagUnmapped(rec) )
		{
			std::string contig = toCString(getContigName(rec, bamFileIn));
			int32_t begin = rec.beginPos;
			int32_t end = rec.beginPos + length(rec.seq);

			String<MyCargo> cargos;
			findIntervals(cargos, trees[contig], begin, end);

			for (int i=0; i<length(cargos); ++i)
			{
				types.insert( cargos[i][0] );
				symbols.insert( cargos[i][2] );
				ids.insert( cargos[i][3] );
			}

			n_genes = symbols.size();

			// NO GENE SO INTERGENIC
			if ( 0 == symbols.size() )
			{
				status = "NO_GENE";
				symbol = "NULL";
				id = "NULL";
				function = "INTERGENIC";
			}

			// ONLY ONE GENE
			else if ( 1 == symbols.size() )
			{
				status = "ONE_GENE";
				it = symbols.begin();
				symbol = toCString(*it);
				it = ids.begin();
				id = toCString(*it);

				int flag = functionFlag(types, FLAGS);

				if ( flag == FUNCTION::GENIC )
				{
					function = "INTRONIC";
				}
				else if ( flag == (FUNCTION::GENIC | FUNCTION::CODING) )
				{
					function = "EXONIC";
				}
				else if ( flag ==
					(FUNCTION::GENIC | FUNCTION::CODING | FUNCTION::UTR) )
				{
					function = "UTR";
				}
				else
				{
					function = "UNKNOWN";
				}

			}

			// several genes
			else
			{
				status = "MULTIPLE_GENES";
				symbol = "NULL";
				id = "NULL";
				function = "NULL";
			}

			types.clear();
			symbols.clear();
			ids.clear();
		}

		// unmapped
		else
		{
			status = "NULL";
			n_genes = 0;
			symbol = "NULL";
			id = "NULL";
			function = "NULL";
		}

		appendTagValue(tagsDict, "gs", status.c_str());
		appendTagValue(tagsDict, "ga", n_genes);
		appendTagValue(tagsDict, "gi", id);
		appendTagValue(tagsDict, "gn", symbol);
		appendTagValue(tagsDict, "gf", function);
		tagsToBamRecord(rec, tagsDict);

		writeRecord(bamFileOut, rec);

		/////////////////////////////////////////////////////////////////////////

		// counter
		counter++;
		if ( counter % (unsigned long)1000000 == 0 ) {
			std::cerr << counter << std::endl;
		}
		//if ( counter > 1000000 ) { break; }
	}

	close(bamFileIn);
	close(bamFileOut);

	////////////////////////////////////////////////////////////////////////////
	return 0;
}
