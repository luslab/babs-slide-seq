#include <map>
#include <set>
#include <seqan/bam_io.h>

using namespace seqan;

typedef std::set<CharString> Umis;
typedef std::map<CharString, Umis> Beads;

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	int threshold = std::stoi(argv[3]);

	unsigned long counter = 0;
	unsigned int tag_idx = 0;

	CharString up_status;
	CharString align_status;
	CharString barcode;
	CharString umi;
	CharString umi_status;

	Beads beads;
	std::map<CharString, CharString> status;

	BamHeader header;
	BamAlignmentRecord rec;

	////////////////////////////////////////////////////////////////////////////

	BamFileIn bam(argv[1]);
	readHeader(header, bam);

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict tagsDict(rec.tags);

		// UP primer status
		if ( findTagKey(tag_idx, tagsDict, "us") ) {
			extractTagValue(up_status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find UP primer status" << std::endl;
			return 1;
		}

		// alignment status
		if ( findTagKey(tag_idx, tagsDict, "as") ) {
			extractTagValue(align_status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find alignment status" << std::endl;
			return 1;
		}

		if ( "MATCHED" == up_status && "MAPPED" == align_status )
		{
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

			if ( beads.find(barcode) == beads.end() )
			{
				Umis umis;
				umis.insert(umi);
				beads[barcode] = umis;
			}
			else
			{
				beads[barcode].insert(umi);
			}
		}

		counter++;
		if ( counter % (unsigned long)1000000 == 0 ) {
			std::cerr << counter << std::endl;
		}
		//if ( counter > 1000000 ) { break; }

	} // bam record

	close(bam);

	////////////////////////////////////////////////////////////////////////////
	// COUNT
	
	std::map<CharString, unsigned long> counts;
	Beads::iterator it;
	for ( it = beads.begin(); it != beads.end(); ++it )
	{
		counts[ it->first ] = it->second.size();
		if ( it->second.size() < threshold )
		{
			status[ it->first ] = "TOO_LOW";
		}
		else
		{
			status[ it->first ] = "PASS";
		}
	}

	////////////////////////////////////////////////////////////////////////////
	// ADD TAG

	BamFileIn bamFileIn(argv[1]);
	BamFileOut bamFileOut(context(bamFileIn), argv[2]);

	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	std::string n_umis;

	counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(rec, bamFileIn);
		BamTagsDict tagsDict(rec.tags);

		// UP primer status
		if ( findTagKey(tag_idx, tagsDict, "us") ) {
			extractTagValue(up_status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find UP primer status" << std::endl;
			return 1;
		}

		// alignment status
		if ( findTagKey(tag_idx, tagsDict, "as") ) {
			extractTagValue(align_status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find alignment status" << std::endl;
			return 1;
		}

		if ( "MATCHED" == up_status && "MAPPED" == align_status )
		{
			// barcode 
			if ( findTagKey(tag_idx, tagsDict, "bc") ) {
				extractTagValue(barcode, tagsDict, tag_idx);
			}
			else
			{
				std::cerr << "Problem: cannot find barcode" << std::endl;
				return 1;
			}

			n_umis = std::to_string( counts[barcode] );
			umi_status = status[barcode];
		}
		else
		{
			n_umis = "NULL";
			umi_status = "NULL";
		}

		appendTagValue(tagsDict, "bn", n_umis);
		appendTagValue(tagsDict, "bt", umi_status);
		tagsToBamRecord(rec, tagsDict);
		writeRecord(bamFileOut, rec);

		/////////////////////////////////////////////////////////////////////////
		counter++;
		if ( counter % (unsigned long)1000000 == 0 ) {
			std::cerr << counter << std::endl;
		}
		//if ( counter > 1000000 ) { break; }

	} // bam record

	close(bamFileIn);
	close(bamFileOut);

	//////////////////////////////////////////////////////////////////////////////
	return 0;
}
