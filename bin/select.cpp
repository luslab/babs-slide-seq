#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <seqan/basic.h>
#include <seqan/bam_io.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////
typedef std::tuple<CharString, CharString> Molecule;
typedef
std::tuple<
	unsigned long, // position in bam file
	CharString, // read name
	CharString, // gene id
	CharString, // function
	int // alignment score
> Mapping;
typedef std::tuple<CharString, int> GeneMapping;
typedef std::tuple<unsigned long, int> FunctionMapping;
typedef std::map<Molecule, std::set<Mapping>> Mappings;

///////////////////////////////////////////
bool compareMapping(Mapping m1, Mapping m2)
{
	int score1 = std::get<4>(m1);
	int score2 = std::get<4>(m2);
	return (score1>score2);
}

/////////////////////////////////////////////////////////
bool compareGeneMapping(GeneMapping gm1, GeneMapping gm2)
{
	int score1 = std::get<1>(gm1);
	int score2 = std::get<1>(gm2);
	return (score1>score2);
}

/////////////////////////////////////////////////////////////////////
bool compareFunctionMapping(FunctionMapping fm1, FunctionMapping fm2)
{
	int rank1 = std::get<1>(fm1);
	int rank2 = std::get<1>(fm2);
	return (rank1<rank2);
}

///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
std::set<GeneMapping> getGeneMappings(std::set<Mapping> mappings)
{
	std::set<Mapping>::iterator it;
	std::set<GeneMapping> gms;

	for ( it = mappings.begin(); it != mappings.end(); ++it )
	{
		GeneMapping gm = std::make_tuple(std::get<2>(*it), std::get<4>(*it));
		gms.insert(gm);
	}

	return gms;
}

///////////////////////////////////////////
int getMaxScore(std::set<Mapping> mappings)
{
	std::set<GeneMapping> gms = getGeneMappings(mappings);
	std::vector<GeneMapping> v(gms.begin(), gms.end());
	std::sort(v.begin(), v.end(), compareGeneMapping);

	int max_score = std::get<1>(v[0]);

	int n = 0;

	for (int i=0; i<v.size(); i++)
	{
		if ( max_score == std::get<1>(v[i]) )
		{
			n++;
		}
	}

	if ( n == 1 )
	{
		return max_score;
	}
	else
	{
		return -1;
	}
}

//////////////////////////////
std::tuple<unsigned long, int>
getMaxScoreRead(std::set<Mapping> mappings, int max_score)
{
	std::set<Mapping>::iterator it;

	unsigned long position;
	int n = 0;

	for ( it = mappings.begin(); it != mappings.end(); ++it )
	{
		int score = std::get<4>(*it);

		if ( max_score == score )
		{
			n++;
			position = std::get<0>(*it);
		}
	}

	return std::make_tuple(position, n);
}

///////////////////////////////////////////////////////////////////////
std::set<Mapping> selectMappings(std::set<Mapping> mappings, int score)
{
	std::set<Mapping> s;
	std::set<Mapping>::iterator it;
	for ( it = mappings.begin(); it != mappings.end(); ++it )
	{
		if ( score == std::get<4>(*it) )
		{
			s.insert(*it);
		}
	}
	return s;
}

////////////////////////////////////////////////////////////////
unsigned long chooseReadFromFunction(std::set<Mapping> mappings)
{
	std::vector<FunctionMapping> v;
	std::set<Mapping>::iterator it;

	std::map<CharString, int> ranks;
	ranks["CODING"] = 0;
	ranks["UTR"] = 1;
	ranks["INTRONIC"] = 2;

	for ( it = mappings.begin(); it != mappings.end(); ++it )
	{
		int rank = ranks[ std::get<3>(*it) ];
		v.push_back( std::make_tuple(std::get<0>(*it), rank) );
	}

	// mappings are already orered by genome position because the bam file is
	// supposed to be sorted because of Picar MarkDuplicates
	std::sort(v.begin(), v.end(), compareFunctionMapping);

	return std::get<0>(v[0]); 
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	unsigned long counter = 0;

	BamFileIn bam(argv[1]);
	BamHeader header;
	readHeader(header, bam);

	BamAlignmentRecord rec;

	CharString read;
	CharString barcode;
	CharString umi;
	CharString multimap;
	CharString id;
	CharString function;
	int score;

	unsigned int tag_idx = 0;

	Mappings mappings;

	std::map<unsigned long, CharString> tag;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict tagsDict(rec.tags);

		read = rec.qName;

		// barcode 
		if ( findTagKey(tag_idx, tagsDict, "bc") ) {
			extractTagValue(barcode, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find barcode" << std::endl;
			return 1;
		}

		// umi
		if ( findTagKey(tag_idx, tagsDict, "mi") ) {
			extractTagValue(umi, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find UMI" << std::endl;
			return 1;
		}

		// multimapping status
		if ( findTagKey(tag_idx, tagsDict, "mm") ) {
			extractTagValue(multimap, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find multi-mapping status" << std::endl;
			return 1;
		}

		// gene id
		if ( findTagKey(tag_idx, tagsDict, "qi") ) {
			extractTagValue(id, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find gene ID" << std::endl;
			return 1;
		}

		// function
		if ( findTagKey(tag_idx, tagsDict, "qf") ) {
			extractTagValue(function, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find function" << std::endl;
			return 1;
		}

		// score 
		if ( findTagKey(tag_idx, tagsDict, "AS") ) {
			extractTagValue(score, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find score" << std::endl;
			return 1;
		}

		Molecule mol  = std::make_tuple(barcode, umi);
		Mapping mapping = std::make_tuple(counter, read, id, function, score);

		// we take only those reads
		if ( "UNIQUE" == multimap || "INCLUDED" == multimap  )
		{
			if ( mappings.find(mol) == mappings.end() )
			{
				std::set<Mapping> s;
				s.insert(mapping);
				mappings[mol] = s;
			}
			else
			{
				mappings[mol].insert(mapping);
			}
		}

		// we ignore ambiguous reads and duplicates
		else
		{
			tag[counter] = "IGNORED";
		}

		/////////////////////////////////////////////////////////////////////////
		// counter
		counter++;
		if ( counter % (unsigned long)1000000 == 0 ) {
			std::cerr << counter << std::endl;
		}
		//if ( counter > 1000000 ) { break; }

	} // bam record

	close(bam);

	////////////////////////////////////////////////////////////////////////////
	
	Mappings::iterator it;

	unsigned long n_unique = 0;
	unsigned long n_chosen = 0;
	unsigned long n_unresolved = 0;

	for ( it = mappings.begin(); it != mappings.end(); ++it )
	{
		std::set<Mapping> s = it->second;
		std::vector<Mapping> v(s.begin(), s.end());

		unsigned long position; // in bam file
		bool determined;

		// only one read: good!
		if ( v.size() == 1 )
		{
			position = std::get<0>(v[0]);
			determined = true;
		}

		// multiple reads
		else
		{
			std::sort(v.begin(), v.end(), compareMapping);
			int max_score = getMaxScore(s);

			// we can't determine unfortunately
			if ( max_score == -1 )
			{
				determined = false;
			}

			// we found a gene that stands out
			else
			{
				determined = true;
				std::tuple<unsigned long, int> msr = getMaxScoreRead(s, max_score);

				// we have only one read
				if ( std::get<1>(msr) == 1 )
				{
					position = std::get<0>(msr);
				}

				// we still have multiple reads, so select by function and then
				// by mapping position on the genome because the bam file is
				// supposed to be sorted already because we called Picard
				// MarkDuplicates previously
				else
				{
					std::set<Mapping> good_mappings  = selectMappings(s, max_score);
					position = chooseReadFromFunction(good_mappings);
				}
			} // end selection by function
		} // multiple reads

		///////////////////////////////////////////////////////////////////////
		// NOW WE ADD THE TAGS

		if ( v.size() == 1 )
		{
			tag[ std::get<0>(v[0]) ] = "UNIQUE";
			n_unique++;
		}
		else
		{
			if ( determined )
			{
				for (int i=0; i<v.size(); i++)
				{
					unsigned long pos = std::get<0>(v[i]);
					
					if ( pos == position )
					{
						tag[pos] = "INCLUDED";
					}
					else
					{
						tag[pos] = "EXCLUDED";
					}
				}
				n_chosen++;
			}
			else
			{
				for (int i=0; i<v.size(); i++)
				{
					tag[ std::get<0>(v[i]) ] = "UNRESOLVED";
				}
				n_unresolved++;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	BamFileIn bamFileIn(argv[1]);
	BamFileOut bamFileOut(context(bamFileIn), argv[2]);
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(rec, bamFileIn);
		BamTagsDict tagsDict(rec.tags);
		appendTagValue(tagsDict, "cs", tag[counter]);
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
