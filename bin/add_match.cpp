#include <map>
#include <set>
#include <seqan/bam_io.h>

using namespace seqan;

int main(int argc, char **argv)
{
	////////////////////////////////////////////////////////////////////////////
	// BARCODE MAPPING

	std::map<CharString, CharString> map;

	std::ifstream csv(argv[1]);
	std::string line;
	CharString seq_barcode, puck_barcode;
	std::set<int> lengths;
	std::set<int>::iterator lit;
	int i = 0;


	while( std::getline(csv, line) )
	{
		// header
		if ( i == 0 )
		{
			i++;
			continue;
		}

		seq_barcode = line.substr(0, line.find(","));
		puck_barcode = line.substr(line.find(",")+1, line.length());

		if ( length(seq_barcode) != length(puck_barcode) )
		{
			std::cerr
				<< "Error: the barcodes pair " << i << " have different lengths"
				<< std::endl;
			return 1;
		}
		else
		{
			lengths.insert( length(seq_barcode) );
			map[seq_barcode] = puck_barcode;
		}
	}
	csv.close();


	int length;

	if ( 1 != lengths.size() )
	{
		std::cerr
			<< "Error: not all the barcodes pair have the same length"
			<< std::endl;
	}
	else
	{
		lit = lengths.begin();
		length = *lit;
	}

	CharString dummy = "";
	for (int i=0; i<length; i++)
	{
		dummy += 'X';
	}

	////////////////////////////////////////////////////////////////////////////
	// ADD TAG

	BamFileIn bamFileIn(argv[2]);
	BamFileOut bamFileOut(context(bamFileIn), argv[3]);

	BamHeader header;
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	BamAlignmentRecord record;
	unsigned int tag_idx = 0;
	CharString barcode;
	CharString bead;
	CharString up_status;
	std::string status;

	long unsigned int counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(record, bamFileIn);

		BamTagsDict tagsDict(record.tags);

		// function
		if ( findTagKey(tag_idx, tagsDict, "bc") )
		{
			extractTagValue(barcode, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find barcode" << std::endl;
			return 1;
		}

		if ( findTagKey(tag_idx, tagsDict, "us") )
		{
			extractTagValue(up_status, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find UP primer status" << std::endl;
			return 1;
		}

		if ( up_status == "MATCHED" )
		{
			if ( map.find(barcode) != map.end() )
			{
				status = "MATCHED";
				bead = map[barcode];
			}
			else
			{
				status = "UNMATCHED";
				bead = dummy;
			}
		}
		else
		{
				status = "NULL";
				bead = dummy;
		}

		appendTagValue(tagsDict, "bs", status.c_str());
		appendTagValue(tagsDict, "bm", toCString(bead));
		tagsToBamRecord(record, tagsDict);

		/////////////////////////////////////////////////////////////////////////

		writeRecord(bamFileOut, record);

		/////////////////////////////////////////////////////////////////////////

		counter++;

		if ( counter % (long unsigned int)1000000 == 0 )
		{
			std::cout << counter << std::endl;
		}
	}

	close(bamFileIn);
	close(bamFileOut);

	////////////////////////////////////////////////////////////////////////////

	return 0;
}
