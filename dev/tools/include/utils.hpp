/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <string>
#include <map>
#include <fstream>

#include <seqan/sequence.h>
#include <seqan/bam_io.h>

// ============================================================================
// Prototypes
// ============================================================================

int HammingDist(std::string, std::string);
std::string DummySeq(int);
int GetPosition(seqan::StringSet<seqan::CharString>, seqan::CharString);

// ----------------------------------------------------------------------------
// WriteCounter()
// ----------------------------------------------------------------------------

template<typename T, typename U>
void WriteCounter(std::string sample, std::string path, std::map<T, U> counter)
{
	std::ofstream csv(path);

	typename std::map<T, U>::iterator it;
	
	for (it = counter.begin(); it != counter.end(); ++it)
	{
		csv << sample << "," << it->first << "," << it->second << std::endl;
	}

	csv.close();
}

// ----------------------------------------------------------------------------
// ExtractTag()
// ----------------------------------------------------------------------------

template<typename T>
void ExtractTag(char* tag, seqan::BamTagsDict dict, T& value)
{
	try
	{
		unsigned int idx = 0;

		if ( seqan::findTagKey(idx, dict, tag) )
		{
			seqan::extractTagValue(value, dict, idx);
		}
		else
		{
			throw(tag);
		}
	}

	catch(char* name)
	{
		std::cerr << "Problem: cannot find " << std::string(name) << " tag" << std::endl;
		exit(EXIT_FAILURE);
	}
}

#endif

