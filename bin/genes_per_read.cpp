#include <map>
#include <set>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>

using namespace seqan;


///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	unsigned long counter = 0;

	BamFileIn bam(argv[1]);
	BamHeader header;
	readHeader(header, bam);

	BamAlignmentRecord rec;

	CharString read;
	CharString id;

	unsigned int tag_idx = 0;

	std::map<CharString, std::set<CharString>> genes;

	std::map<unsigned long, CharString> tag;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict tagsDict(rec.tags);

		read = rec.qName;

		// gene ID
		if ( findTagKey(tag_idx, tagsDict, "qi") ) {
			extractTagValue(id, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find gene ID" << std::endl;
			return 1;
		}

		if ( genes.find(read) == genes.end() )
		{
			std::set<CharString> s;
			s.insert(id);
			genes[read] = s;
		}
		else
		{
			genes[read].insert(id);
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
	////////////////////////////////////////////////////////////////////////////
	
	std::map<CharString, std::set<CharString>>::iterator it;

	std::map<CharString, CharString> status;

	unsigned long unassigned = 0;
	unsigned long one_gene = 0;
	unsigned several_genes = 0;

	std::set<CharString> null_null_set;
	null_null_set.insert("NULL");

	for ( it = genes.begin(); it != genes.end(); ++it )
	{
		if ( it->second == null_null_set )
		{
			status[ it->first ] = "UNASSIGNED";
			unassigned++;
		}

		if ( it->second != null_null_set && it->second.size() == 1 )
		{
			status[ it->first ] = "ONE_GENE";
			one_gene++;

		}

		if ( it->second.size() > 1 )
		{
			status[ it->first ] = "SEVERAL_GENES";
			several_genes++;

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
		appendTagValue(tagsDict, "rm", status[rec.qName]);
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

