#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <seqan/basic.h>
#include <seqan/bam_io.h>

using namespace seqan;

///////////////////////////////////////////////////////
bool contains(std::set<CharString> tags, CharString tag)
{
	if ( tags.find(tag) == tags.end() )
	{
		return false;
	}
	else
	{
		return true;
	}
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	unsigned long counter = 0;

	BamFileIn bam(argv[1]);
	BamHeader header;
	readHeader(header, bam);
	BamAlignmentRecord rec;

	unsigned int tag_idx = 0;

	CharString read;
	CharString status;
	std::map<CharString, std::set<CharString>> tags;
	std::map<CharString, std::set<CharString>>::iterator it;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict tagsDict(rec.tags);

		read = rec.qName;

		// barcode 
		if ( findTagKey(tag_idx, tagsDict, "dp") ) {
			extractTagValue(status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find duplicate status" << std::endl;
			return 1;
		}

		if ( tags.find(read) == tags.end() )
		{
			std::set<CharString> s;
			s.insert(status);
			tags[read] = s;
		}
		else
		{
			tags[read].insert(status);
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
	
	std::map<CharString, CharString> annot;
	
	for ( it = tags.begin(); it != tags.end(); ++it )
	{
		CharString read = it->first;
		std::set<CharString> s = it->second;
		bool is_primary = contains(s, "PRIMARY");
		bool is_duplicate = contains(s, "DUPLICATE");
		bool is_null = contains(s, "NULL");

		if ( is_primary && !is_null )
		{
			annot[read] = "PRIMARY";
		}
		else if ( !is_primary && is_duplicate && !is_null )
		{
			annot[read] = "DUPLICATE";
		}
		else if ( !is_primary && !is_duplicate && is_null )
		{
			annot[read] = "NULL";
		}
		else
		{
			std::cerr << "Error: read with impossible combination" << std::endl;
			return 1;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////

	BamFileIn bamFileIn(argv[1]);
	BamFileOut bamFileOut(context(bamFileIn), argv[2]);
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(rec, bamFileIn);
		read = rec.qName;
		BamTagsDict tagsDict(rec.tags);
		appendTagValue(tagsDict, "ds", annot[read]);
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

