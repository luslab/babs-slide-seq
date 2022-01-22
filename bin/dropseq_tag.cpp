#include <map>
#include <set>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>

using namespace seqan;

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
int main(int argc, char **argv)
{
	////////////////////////////////////////////////////////////////////////////

	GffFileIn gff(argv[1]);
	GffRecord record;
	
	int g = 1;

	std::map<CharString, std::tuple<CharString, CharString>> transcripts;

	while ( ! atEnd(gff) )
	{
		readRecord(record, gff);

		int namePos = get_position(record.tagNames, "gene_name");
		int idPos = get_position(record.tagNames, "gene_id");
		int transPos = get_position(record.tagNames, "transcript_id");

		if ( -1 != namePos && -1 != idPos && -1 != transPos )
		{
			std::tuple<CharString, CharString> gene = 
				std::make_tuple(
					record.tagValues[idPos],
					record.tagValues[namePos]
			);

			transcripts[ record.tagValues[transPos] ] = gene;

		}
		// counter
		if ( g % (unsigned long)100000 == 0 ) {
			std::cerr << g << " gff entries" << std::endl;
		}
		g++;
	}

	close(gff);

	////////////////////////////////////////////////////////////////////////////

	BamFileIn bamFileIn(argv[2]);
	BamFileOut bamFileOut(context(bamFileIn), argv[3]);

	BamHeader header;
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	BamAlignmentRecord rec;

	std::tuple<CharString, CharString> gene;
	CharString function;
	CharString transcript;
	CharString id;
	CharString name;

	unsigned long counter = 0;
	unsigned int tag_idx = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(rec, bamFileIn);
		BamTagsDict tagsDict(rec.tags);

		if ( findTagKey(tag_idx, tagsDict, "XF") ) {
			extractTagValue(function, tagsDict, tag_idx);
		}
		else
		{
			function = "NULL";
		}

		if ( findTagKey(tag_idx, tagsDict, "GE") ) {
			extractTagValue(transcript, tagsDict, tag_idx);
			gene = transcripts[transcript];
			id = std::get<0>(gene);
			name = std::get<1>(gene);
		}
		else
		{
			transcript = "NULL";
			id = "NULL";
			name = "NULL";
		}

		appendTagValue(tagsDict, "qi", id);
		appendTagValue(tagsDict, "qn", name);
		appendTagValue(tagsDict, "qf", function);
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
