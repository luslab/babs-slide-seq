/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>

#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "gene.hpp"
#include "counts.hpp"
#include "counter.hpp"

int main(int argc, char** argv)
{

	////////////////////////////////////////////////////////////////////////////
	std::set<Gene> genes;
	std::vector<Gene> genes_v;

	for (int i=0; i<26; i++)
	{
		char* id;
		char* symbol;
		id = (char*)malloc(sizeof(char)*9);
		symbol = (char*)malloc(sizeof(char)*6);

		sprintf(id, "ENSMUS%0.2d", 1+i);
		sprintf(symbol, "Gene%c", 90-i);

		genes.insert( Gene(id, symbol) );
		genes_v.push_back( Gene(id, symbol) );
	}

	std::cout << std::endl;
	std::cout << "The genes:" << std::endl;
	for (auto& gene : genes)
	{
		std::cout << gene << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////
	
	srand(123);
	
	char BASE[4] = {'A', 'C', 'G', 'T'};

	std::set<std::string> barcodes;
	std::vector<std::string> barcodes_v;

	for (int i=0; i<4; i++)
	{
		char sequence[14];
		for (int j=0; j<14; j++)
		{
			sequence[j] = BASE[ rand() % 4 ];
		}

		
		barcodes.insert( std::string(sequence) );
		barcodes_v.push_back( std::string(sequence) );
	}

	std::cout << std::endl;
	std::cout << "The barcodes:" << std::endl;
	for (auto& barcode : barcodes)
	{
		std::cout << barcode << std::endl;
	}
	
	////////////////////////////////////////////////////////////////////////////

	Counter counter = Counter(genes);

	for (int i=0; i<26; i++)
	{
		char* id = (char*)genes_v[i].GetID().c_str();
		if ( (i+1) % 2 == 0 )
		{
			counter.Increment((char*)barcodes_v[0].c_str(), id);
		}

		else if ( (i+1) % 3 == 0 )
		{
			counter.Increment((char*)barcodes_v[0].c_str(), id);
			counter.Increment((char*)barcodes_v[0].c_str(), id);
			counter.Increment((char*)barcodes_v[1].c_str(), id);
		}

		else if ( (i+1) % 5 == 0 )
		{
			counter.Increment((char*)barcodes_v[2].c_str(), id);
			counter.Increment((char*)barcodes_v[2].c_str(), id);
			counter.Increment((char*)barcodes_v[2].c_str(), id);
		}
	}

	std::cout << std::endl;
	std::cout << "Test of the Increment() method:" << std::endl;
	std::cout << counter << std::endl;

	std::cout << std::endl;
	std::cout << "Test of the GetNumberOf*() methods:" << std::endl;
	std::cout << "genes: " << counter.GetNumberOfGenes() << std::endl;
	std::cout << "barcodes: " << counter.GetNumberOfBarcodes() << std::endl;
	std::cout << "entries: " << counter.GetNumberOfEntries() << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the GetGenes() method:" << std::endl;
	for (auto& gene : counter.GetGenes())
	{
		std::cout << gene << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Test of the GetBarcodes() method:" << std::endl;
	for (auto& barcode : counter.GetBarcodes())
	{
		std::cout << barcode << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the GetDGE() method:" << std::endl;
	counter.GetDGE(std::cout, std::cout, std::cout);

	std::ofstream feat("features.tsv");
	std::ofstream bcd("barcodes.tsv");
	std::ofstream mtx("matrix.mtx");
	counter.GetDGE(feat, bcd, mtx);
	feat.close();
	bcd.close();
	mtx.close();

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

