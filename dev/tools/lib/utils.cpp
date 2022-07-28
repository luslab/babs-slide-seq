/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <map>

#include <seqan/sequence.h>

// ----------------------------------------------------------------------------
// HammingDist()
// ----------------------------------------------------------------------------

int HammingDist(std::string sequence1, std::string sequence2)
{
	int distance = 0;

	for (int i=0; i<sequence1.size(); i++)
	{
		char base1 = sequence1[i];
		char base2 = sequence2[i];

		if ( base1 != base2 )
		{
			distance++;
		}
	}

	return distance;
}

// ----------------------------------------------------------------------------
// DummySeq()
// ----------------------------------------------------------------------------

std::string DummySeq(int length)
{
	std::string sequence = "";

	for (int i=0; i<length; i++)
	{
		sequence += 'X';
	}

	return sequence;
}

// ----------------------------------------------------------------------------
// GetPosition()
// ----------------------------------------------------------------------------

int GetPosition(seqan::StringSet<seqan::CharString> tag_names, seqan::CharString tag_name)
{
	int position = -1;

	for (int i=0; i<seqan::length(tag_names); i++)
	{
		if ( tag_name == tag_names[i] )
		{
			position = i;
		}
	}

	return position;
}

