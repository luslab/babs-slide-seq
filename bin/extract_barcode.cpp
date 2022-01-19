#include <map>
#include <set>
#include <fstream>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

///////////////////////////////////////////////////////////////////////////////

typedef std::map<char, std::vector< std::tuple<int, int> > > Segments;

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////
int hammingDist(std::string seq1, std::string seq2)
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

////////////////////////////////////////////
int get_length(Segments segments, char type)
{
	int length = 0;

	for (int i=0; i<segments[type].size(); i++)
	{
		length += std::get<1>(segments[type][i]) - std::get<0>(segments[type][i]);
	}

	return length;
}

//////////////////////////////////////
std::string dummy_sequence(int length)
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
	SeqFileIn fastq1(argv[1]);
	SeqFileIn fastq2(argv[2]);
	std::string structure(argv[3]);
	SeqFileOut fastq(argv[4]);
	int up_thresh = atoi(argv[5]);
	std::string sample(argv[6]);

	int length_thresh = 41;

	CharString meta1;
	CharString meta2;
	CharString read1;
	CharString read2;
	CharString qual1;
	CharString qual2;

	////////////////////////////////////////////////////////////////////////////

	Segments segments;
	segments['C'] = std::vector< std::tuple<int, int> >();
	segments['U'] = std::vector< std::tuple<int, int> >();
	segments['X'] = std::vector< std::tuple<int, int> >();
	segments['M'] = std::vector< std::tuple<int, int> >();

	char type;
	int begin = 0;
	int end = 0;
	int size = 0;
	std::string size_str = "";

	for (int i=0; i<structure.size(); i++)
	{
		if ( isdigit(structure[i]) )
		{
			size_str += structure[i];
		}
		else
		{
			type = structure[i];
			size = atoi( size_str.c_str() );
			size_str = "";
			end = begin + size;
			segments[type].push_back( std::tuple<int, int>{begin, end} );
			begin = end;
		}
	}

	int up_length = get_length(segments, 'U');
	int barcode_length = get_length(segments, 'C');
	int umi_length = get_length(segments, 'M');

	std::string dummy_up = dummy_sequence(up_length);
	std::string dummy_barcode = dummy_sequence(barcode_length);
	std::string dummy_umi = dummy_sequence(umi_length);

	////////////////////////////////////////////////////////////////////////////

	long unsigned int total = 0;
	long unsigned int too_short = 0;
	long unsigned int long_enough = 0;
	long unsigned int up_match = 0;
	long unsigned int up_nonmatch = 0;

	////////////////////////////////////////////////////////////////////////////
	
	std::map<int, long unsigned int> distances;
	std::map<int, long unsigned int>::iterator it;
	distances[-1] = 0;
	for (int i=0; i<=up_length; i++)
	{
		distances[i] = 0;
	}

	////////////////////////////////////////////////////////////////////////////

	std::string length_status;
	std::string up;
	int up_hamming;
	std::string up_status;
	std::string barcode;
	std::string umi;

	while ( ! atEnd(fastq1) )
	{
		readRecord(meta1, read1, qual1, fastq1);
		readRecord(meta2, read2, qual2, fastq2);

		if ( length(read1) < length_thresh )
		{
			too_short++;

			length_status = "TOO_SHORT";
			up = dummy_up;
			up_hamming = -1;
			up_status = "UNMATCHED";
			barcode = dummy_barcode;
			umi = dummy_umi;
		}
		else
		{
			long_enough++;

			length_status = "LONG_ENOUGH";

			std::map<char, CharString> sequences;
			sequences['C'] = ""; // bead
			sequences['U'] = ""; // UP primer
			sequences['X'] = "";
			sequences['M'] = ""; // UMI
			for (Segments::iterator it=segments.begin(); it!=segments.end(); ++it)
			{
				for (int i=0; i<it->second.size(); i++)
				{
					int begin = std::get<0>(it->second[i]);
					int end = std::get<1>(it->second[i]);
					sequences[it->first] += infix(read1, begin, end);
				}
			}

			up = toCString(sequences['U']);
			barcode = toCString(sequences['C']);
			umi = toCString(sequences['M']);

			up_hamming = hammingDist(up, "TCTTCAGCGTTCCCGAGA");

			if ( up_hamming <= up_thresh )
			{
				up_status = "MATCHED";
				up_match++;
			}
			else
			{
				up_status = "UNMATCHED";
				up_nonmatch++;
			}
		}

		distances[up_hamming]++;

		/////////////////////////////////////////////////////////////////////////
		CharString meta;
		append(meta, length_status);
		append(meta, ":");
		append(meta, up);
		append(meta, ":");
		append(meta, std::to_string(up_hamming));
		append(meta, ":");
		append(meta, up_status);
		append(meta, ":");
		append(meta, barcode);
		append(meta, ":");
		append(meta, umi);
		append(meta, ":");
		append(meta, meta2);
		
		writeRecord(fastq, meta, read2, qual2);

		/////////////////////////////////////////////////////////////////////////

		total++;
		if ( total % (long unsigned int)1000000 == 0 )
		{
			std::cout << total << std::endl;
		}
	}

	close(fastq1);
	close(fastq2);
	close(fastq);

	////////////////////////////////////////////////////////////////////////////
	std::ofstream metrics( sample + ".extract_barcode.csv" );
	metrics << sample << ",UP primer length," << up_length << std::endl;
	metrics << sample << ",Bead barcode length," << barcode_length << std::endl;
	metrics << sample << ",UMI length," << umi_length << std::endl;
	metrics << sample << ",Total reads," << total << std::endl;
	metrics << sample << ",Too short reads," << too_short << std::endl;
	metrics << sample << ",Long enough reads," << long_enough << std::endl;
	metrics << sample << ",UP primer match," << up_match << std::endl;
	metrics << sample << ",UP primer non match," << up_nonmatch << std::endl;
	metrics.close();

	////////////////////////////////////////////////////////////////////////////
	std::ofstream dists( sample + ".extract_barcode.up_distances.csv" );
	for (it=distances.begin(); it!=distances.end(); ++it)
	{
		dists << sample << "," << it->first << "," << it->second << std::endl;
	}
	dists.close();

	////////////////////////////////////////////////////////////////////////////
	return 0;
}
