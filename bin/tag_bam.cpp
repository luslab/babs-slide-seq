#include <regex>
#include <map>
#include <tuple>
#include <algorithm>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////

typedef std::map<char, std::vector< std::tuple<int, int> > > Segments;

///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
int hammingDist(std::string seq1, std::string seq2)//
/////////////////////////////////////////////////////
{
	int distance = 0;

	for (int i=0; i<seq1.size(); i++)
	{
		char s1 = seq1[i];
		char s2 = seq2[i];

		if ( s1 != s2 )
		{
			distance++;
		}
	}

	return distance;
}

//////////////////////////////////////////////
int get_length(Segments segments, char type)//
//////////////////////////////////////////////
{
	int length = 0;

	for (int i=0; i<segments[type].size(); i++)
	{
		length += std::get<1>(segments[type][i]) - std::get<0>(segments[type][i]);
	}

	return length;
}

////////////////////////////////////////
std::string dummy_sequence(int length)//
////////////////////////////////////////
{
	std::string sequence = "";

	for (int i=0; i<length; i++)
	{
		sequence += 'X';
	}

	return sequence;
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	BamFileIn bamFileIn(argv[1]);
	BamFileOut bamFileOut(context(bamFileIn), argv[2]);

	BamHeader header;
	readHeader(header, bamFileIn);
	writeHeader(bamFileOut, header);

	BamAlignmentRecord record;

	////////////////////////////////////////////////////////////////////////////

	std::string regex;
	regex.append("([^:]+)"); // length_status
	regex.append(":");
	regex.append("([A-Z]+)"); // up
	regex.append(":");
	regex.append("(-*[0-9]+)"); // up_hamming
	regex.append(":");
	regex.append("([^:]+)"); // up_status
	regex.append(":");
	regex.append("([A-Z]+)"); // barcode
	regex.append(":");
	regex.append("([A-Z]+)"); // umi
	regex.append(":");
	regex.append("(.*)$"); // read_name
	std::regex rgx(regex);

	std::smatch m;

	////////////////////////////////////////////////////////////////////////////

	std::string length_status;
	std::string up;
	int up_hamming;
	std::string up_status;
	std::string barcode;
	std::string umi;
	std::string read_name;
	std::string mapping_status;
	std::string duplicate_status;

	////////////////////////////////////////////////////////////////////////////
	
	long unsigned int counter = 0;

	while ( ! atEnd(bamFileIn) )
	{
		readRecord(record, bamFileIn);

		/////////////////////////////////////////////////////////////////////////

		std::string qname( toCString(record.qName) );

		if ( std::regex_match(qname, m, rgx) )
		{
			length_status = m[1];
			up = m[2];
			up_hamming = std::stoi(m[3]);
			up_status = m[4];
			barcode = m[5];
			umi = m[6];
			read_name = m[7];
		}
		else
		{
			std::cerr << "Unable to match regex" << std::endl;
			std::cerr << regex << std::endl;
			std::cerr << qname << std::endl;
			return 1;
		}

		/////////////////////////////////////////////////////////////////////////

		record.qName = read_name.c_str();

		if ( hasFlagUnmapped(record) )
		{
			mapping_status = "UNMAPPED";
			duplicate_status = "NULL";
		}
		else
		{
			mapping_status = "MAPPED";

			if ( hasFlagDuplicate(record) )
			{
				duplicate_status = "DUPLICATE";
			}
			else
			{
				duplicate_status = "PRIMARY";
			}
		}

		BamTagsDict tagsDict(record.tags);
		appendTagValue(tagsDict, "ls", length_status.c_str());
		appendTagValue(tagsDict, "up", up.c_str());
		appendTagValue(tagsDict, "uh", up_hamming);
		appendTagValue(tagsDict, "us", up_status.c_str());
		appendTagValue(tagsDict, "bc", barcode.c_str());
		appendTagValue(tagsDict, "mi", umi.c_str());
		appendTagValue(tagsDict, "as", mapping_status.c_str());
		appendTagValue(tagsDict, "dp", duplicate_status.c_str());
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

	return 0;
}
