/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <set>
#include <iostream>

#include "gene.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	
	Gene gene1 = Gene("ENSMUS000", "GeneA");
	Gene gene2 = Gene("ENSMUS001", "GeneB");
	Gene gene3 = Gene("ENSMUS002", "GeneC");
	Gene gene4 = Gene("ENSMUS002", "GeneC");
	Gene gene5 = Gene("ENSMUS002", "GeneD");

	std::set<Gene> genes;
	genes.insert(gene1);
	genes.insert(gene2);
	genes.insert(gene3);
	genes.insert(gene4);
	genes.insert(gene5);

	std::cout << std::endl;
	std::cout << "Test for Gene in set:" << std::endl;
	for (auto& gene : genes)
	{
		std::cout << gene << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}



