/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdlib.h>
#include <time.h>

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <fstream>
#include <regex>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/gff_io.h>

#include "data.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Data::Data()
{
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Data::Data(char* fasta, char* fasta_index, char* gff, char* puck, unsigned long long n)
{
	this->gff = std::string(gff);
	this->fasta = std::string(fasta);
	this->fasta_index = std::string(fasta_index);
	this->features = Features();
	this->features = GetFeatures(gff);

	// Fake barcodes
	std::set<std::string> seq;
	while ( seq.size() < n ) // !! DANGEROUS !!
	{
		seq.insert(RandSeq(15));
	}
	std::vector<std::string> sequences(seq.begin(), seq.end());

	// Fake coordinates
	std::set<std::tuple<double, double>> all_coords = GetCoordsFromFile(puck);
	std::set<std::tuple<double, double>> random_coords = SampleSet(all_coords, n);
	std::vector<std::tuple<double, double>> coords(random_coords.begin(), random_coords.end());

	this->barcodes = std::vector<std::tuple<std::string, double, double>>();
	for (int i=0; i<n; i++)
	{
		this->barcodes.push_back( std::make_tuple(sequences[i], std::get<0>(coords[i]), std::get<1>(coords[i])) );
	}
}

// ============================================================================
// Getter
// ============================================================================

// ----------------------------------------------------------------------------
// GetBarcodes()
// ----------------------------------------------------------------------------

std::vector<std::tuple<std::string, double, double>> Data::GetBarcodes()
{
	return barcodes;
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// RandomSequence()
// ----------------------------------------------------------------------------

std::string Data::RandSeq(int n)
{
	std::string sequence;

	int bases[4] = {'A', 'C', 'G', 'T'};

	for (int i=0; i<n; i++)
	{
		sequence.push_back(bases[ rand() % 4 ]);
	}
	
	return sequence;
}

// ----------------------------------------------------------------------------
// GetCoordsFromFile()
// ----------------------------------------------------------------------------

std::set<std::tuple<double, double>> Data::GetCoordsFromFile(char* path)
{
	std::set<std::tuple<double, double>> coords;
	std::string line;
	std::ifstream file(path);

	std::regex regex("^([0-9.-]+),([0-9.-]+)$");
	std::smatch m;

	while( std::getline(file, line) )
	{
		if ( std::regex_match(line, m, regex) )
		{
			double x = std::stod(m[1]);
			double y = std::stod(m[2]);
			coords.insert( std::make_tuple(x, y) );
		}
	}
	file.close();

	return coords;
}

// ----------------------------------------------------------------------------
// GetFeatures()
// ----------------------------------------------------------------------------

Features Data::GetFeatures(char* gff_path)
{
	seqan::GffFileIn gff(gff_path);
	seqan::GffRecord record;

	std::set<seqan::CharString> TYPES = {
		"gene", "exon",
		"start_codon", "stop_codon",
		"five_prime_utr", "three_prime_utr"
	};

	std::set<seqan::CharString>::iterator it;

	Features feat;

	int counter = 1;

	while ( ! seqan::atEnd(gff) )
	{
		seqan::readRecord(record, gff);

		it = TYPES.find(record.type);
		if ( it == TYPES.end() )
		{
			continue;
		}

		int gene_id_position = -1;
		for (int i=0; i<seqan::length(record.tagNames); i++)
		{
			if ( "gene_id" == record.tagNames[i] )
			{
				gene_id_position = i;
			}
		}
		
		if ( -1 != gene_id_position && record.endPos - record.beginPos > 200 )
		{
			feat.push_back( std::make_tuple(record.ref, record.beginPos, record.endPos) );
		}

		// counter
		if ( counter % (unsigned long)100000 == 0 ) {
			std::cerr << counter << " gff entries" << std::endl;
		}
		counter++;
		if ( counter > 1000 ) { break; }
	}
	
	return feat;
}

// ----------------------------------------------------------------------------
// CreateQual()
// ----------------------------------------------------------------------------

seqan::CharString Data::CreateQual(int n)
{
	seqan::CharString quality;

	for (int i=0; i<n; i++)
	{
		seqan::appendValue(quality, (char) 63 + (rand() % 10));
	}

	return quality;
}

// ----------------------------------------------------------------------------
// AddRand()
// ----------------------------------------------------------------------------

void Data::AddRand(seqan::CharString& sequence, int n)
{
	std::map<char, char> mutations{ {'A', 'C'}, {'C', 'T'}, {'G', 'A'}, {'T', 'G'} };

	std::set<int> positions;
	while ( positions.size() < n )
	{
		int position = rand() % seqan::length(sequence);
		positions.insert(position);
	}

	for (auto& position : positions)
	{
		sequence[position] = mutations[ sequence[position] ];
	}
}

// ----------------------------------------------------------------------------
// AddRand()
// ----------------------------------------------------------------------------

void Data::AddRand(std::string& sequence, int n)
{
	std::map<char, char> mutations{ {'A', 'C'}, {'C', 'T'}, {'G', 'A'}, {'T', 'G'} };

	std::set<int> positions;
	while ( positions.size() < n )
	{
		int position = rand() % seqan::length(sequence);
		positions.insert(position);
	}

	for (auto& position : positions)
	{
		sequence[position] = mutations[ sequence[position] ];
	}
}

// ----------------------------------------------------------------------------
// GetSequences()
// ----------------------------------------------------------------------------

void Data::GetSequences(seqan::SeqFileOut& fastq1, seqan::SeqFileOut& fastq2,
	unsigned long long n, std::string prefix)
{
	seqan::FaiIndex fai;

	if ( ! seqan::open(fai, this->fasta.c_str(), this->fasta_index.c_str()))
	{
		std::cerr << "Error: cannot open the FASTA index" << std::endl;
		exit(EXIT_FAILURE);
	}

	seqan::CharString transcript;
	seqan::CharString read1;
	seqan::CharString read2;

	int counter = 1;

	for (int i=0; i<n; i++)
	{
		/////////////////////////////////////////////////////////////////////////
		// read 1

		std::string up = "TCTTCAGCGTTCCCGAGA";
		AddRand(up, rand() % 6);
		std::tuple<std::string, double, double> bcd = this->barcodes[ rand() % this->barcodes.size() ];
		std::string barcode = std::get<0>(bcd);
		AddRand(barcode, rand() % 2);
		std::string umi = RandSeq(9);
		read1 = barcode.substr(0, 8) + up + barcode.substr(8, 7) + "TC" + umi;
		if ( rand() % 20 == 0 )
		{
			seqan::resize(read1, 0, rand() % seqan::length(read1));
		}

		/////////////////////////////////////////////////////////////////////////
		// read2

		Feature feat = this->features[ rand() % this->features.size() ];
		seqan::CharString ref = std::get<0>(feat);
		unsigned index = 0;
		if ( ! seqan::getIdByName(index, fai, ref) )
		{
			std::cerr << "Error: cannot find for " << ref << std::endl;
			exit(EXIT_FAILURE);
		}

		seqan::readRegion(transcript, fai, index, std::get<1>(feat), std::get<2>(feat));
		seqan::toUpper(transcript);
		seqan::clear(read2);
		for (int j=0; j<50; j++)
		{
			appendValue(read2, transcript[j]);
		}
		AddRand( read2 , rand() % (seqan::length(read2)/3) );

		seqan::writeRecord(fastq1, prefix+"read"+std::to_string(i+1), read1, CreateQual(seqan::length(read1)));
		seqan::writeRecord(fastq2, prefix+"read"+std::to_string(i+1), read2, CreateQual(seqan::length(read2)));

		// counter
		if ( counter % (unsigned long)10000 == 0 ) {
			std::cerr << counter << " fastq entries" << std::endl;
		}
		counter++;
	}
}

