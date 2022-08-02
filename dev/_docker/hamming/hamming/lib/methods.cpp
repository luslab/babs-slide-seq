/*
 * nourdine.bah@crick.ac.uk
 */

#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "types.hpp"

//////////////////////////////////////////////////
std::vector<std::string> get_sequences(char* path)
{
	std::string sequence;
	std::set<std::string> seq;
	std::ifstream file(path);
	while( std::getline(file, sequence) )
	{
		seq.insert(sequence);
	}
	file.close();

	std::vector<std::string> sequences(seq.begin(), seq.end());
	std::sort(sequences.begin(), sequences.end());

	return sequences;
};

/////////////////////////////////////////
Sequence seq_to_num(std::string sequence)
{
	Sequence number = 0;

	for (int i=0; i<sequence.size(); i++)
	{
		int pos = sequence.size() - i - 1;

		Sequence base;

		if ( 'A' == sequence[pos] ) { base = 1; }
		else if ( 'C' == sequence[pos] ) { base = 2; }
		else if ( 'G' == sequence[pos] ) { base = 3; }
		else if ( 'T' == sequence[pos] ) { base = 4; }
		else if ( 'N' == sequence[pos] ) { base = 5; }
		else { base = 0; }

		base = base << (i * 3);
		number += base;
	}

	return number;
};

///////////////////////////////////////////////////
std::string num_to_seq(Sequence number, int length)
{
	Sequence base;
	Sequence mask = 7;
	std::string sequence;
	sequence.resize(length);

	for (int i=0; i<length; i++)
	{
		int pos = length - i - 1;

		base = number & ( mask << ( i * 3 ) );
		base = base >> ( i * 3 );

		if ( 1 == base ) { sequence[pos] = 'A'; }
		else if ( 2 == base ) { sequence[pos] = 'C'; }
		else if ( 3 == base ) { sequence[pos] = 'G'; }
		else if ( 4 == base ) { sequence[pos] = 'T'; }
		else if ( 5 == base ) { sequence[pos] = 'N'; }
		else { sequence[pos] = 'X'; }
	}

	return sequence;
};

//////////////////////////////////////////////////////////////
char distance(Sequence number1, Sequence number2, char length)
{
	Sequence mask = 7;
	char distance = 0;
	Sequence base1;
	Sequence base2;

	for (int i=0; i<length; i++)
	{
		base1 = number1 & ( mask << ( i * 3 ) );
		base2 = number2 & ( mask << ( i * 3 ) );

		if ( base1 != base2 )
		{
			distance++;
		}
	}

	return distance;
};

///////////////
char get_length(
		std::vector<std::string> reads,
		std::vector<std::string> barcodes,
		char* length)
{
	// iterator
	std::set<char>::iterator it;

	// return -1 if multiple read lengths
	std::set<char> read_lengths;
	for (int i=0; i<reads.size(); i++)
	{
		read_lengths.insert( reads[i].size() );
	}
	if ( 1 != read_lengths.size() )
	{
		return -1;
	}
	it = read_lengths.begin();
	char read_length = *it;

	// return -2 if multiple barcode lengths
	std::set<char> barcode_lengths;
	for (int i=0; i<barcodes.size(); i++)
	{
		barcode_lengths.insert( barcodes[i].size() );
	}
	if ( 1 != barcode_lengths.size() )
	{
		return -2;
	}
	it = barcode_lengths.begin();
	char barcode_length = *it;

	// return -3 if lengths are different
	if ( read_length != barcode_length )
	{
		return -1;
	}
	else
	{
		*length = read_length;
		return 0;
	}
}
