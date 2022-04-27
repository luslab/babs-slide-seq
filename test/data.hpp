/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __DATA_HPP
#define __DATA_HPP

#include <stdlib.h>

#include <tuple>
#include <vector>
#include <set>
#include <map>
#include <regex>

#include <seqan/seq_io.h>
#include <seqan/basic.h>

typedef std::tuple<seqan::CharString, int32_t, int32_t> Feature;
typedef std::vector<Feature> Features;

class Data
{
	private:
		std::string fasta;
		std::string fasta_index;
		std::string gff;
		std::vector<std::tuple<std::string, double, double>> barcodes;
		Features features;
		char bases[4];
	
	public:
		Data();
		Data(char*, char*, char*, char*, unsigned long long);
		Features GetFeatures(char*);
		std::vector<std::tuple<std::string, double, double>> GetBarcodes();
		void AddRand(seqan::CharString&, int);
		void AddRand(std::string&, int);
		seqan::CharString CreateQual(int);
		void GetSequences(seqan::SeqFileOut&, seqan::SeqFileOut&, unsigned long long, std::string prefix);
		std::string RandSeq(int);
		std::set<std::tuple<double, double>> GetCoordsFromFile(char*);
};

// ----------------------------------------------------------------------------
// SampleSet()
// ----------------------------------------------------------------------------

template<typename T>
std::set<T> SampleSet(std::set<T> set, unsigned long long n)
{
	if ( n > set.size() )
	{
		std::cerr << "Error the number of elements to sample is higher than the set size" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<T> vector(set.begin(), set.end());
	std::set<T> sampled;

	while ( sampled.size() < n )
	{
		sampled.insert( vector[ rand() % set.size() ] );
	}

	return sampled;
}

#endif

