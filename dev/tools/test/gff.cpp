/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>

#include <set>
#include <string>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/gff_io.h>

#include "gene.hpp"

int main(int argc, char** argv)
{

	////////////////////////////////////////////////////////////////////////////
	std::set<Gene> genes;

	std::string gff_path("/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf");
	seqan::GffFileIn gff(gff_path.c_str());

	seqan::GffRecord record;

	std::set<seqan::CharString> TYPES = {
		"gene", "exon",
		"start_codon", "stop_codon",
		"five_prime_utr", "three_prime_utr"
	};

	std::set<seqan::CharString>::iterator it;

	int counter = 1;

	while ( ! seqan::atEnd(gff) )
	{
		seqan::readRecord(record, gff);

		it = TYPES.find(record.type);
		if ( it == TYPES.end() )
		{
			continue;
		}

		int gene_name_position = -1;
		for (int i=0; i<seqan::length(record.tagNames); i++)
		{
			if ( "gene_name" == record.tagNames[i] )
			{
				gene_name_position = i;
			}
		}

		int gene_id_position = -1;
		for (int i=0; i<seqan::length(record.tagNames); i++)
		{
			if ( "gene_id" == record.tagNames[i] )
			{
				gene_id_position = i;
			}
		}

		if ( -1 != gene_name_position && -1 != gene_id_position )
		{

			Gene gene = Gene(seqan::toCString(record.tagValues[gene_id_position]), seqan::toCString(record.tagValues[gene_name_position]));
			genes.insert(gene);
		}

		// counter
		if ( counter % (unsigned long)100000 == 0 ) {
			std::cerr << counter << " gff entries" << std::endl;
		}
		counter++;
	}

	for (auto& gene : genes)
	{
		std::cout << gene << std::endl;
	}
	
	close(gff);

	////////////////////////////////////////////////////////////////////////////
	return 0;
}

