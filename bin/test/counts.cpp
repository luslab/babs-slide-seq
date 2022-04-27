/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>

#include "gene.hpp"
#include "counts.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	
	Gene gene1 = Gene("ENSMUS000", "GeneA");
	Gene gene2 = Gene("ENSMUS001", "GeneB");
	Gene gene3 = Gene("ENSMUS002", "GeneC");
	Gene gene4 = Gene("ENSMUS002", "GeneC");
	Gene gene5 = Gene("ENSMUS002", "GeneD");
	Gene gene6 = Gene("ENSMUS003", "GeneE");

	Counts counts = Counts(gene1);
	counts.Increment(gene2);
	counts.Increment(gene2);
	counts.Increment(gene2);
	counts.Increment(gene3);

	std::cout << std::endl;
	std::cout << "Test for the GetSize() method:" << std::endl;
	std::cout << counts.GetSize() << std::endl;

	std::cout << std::endl;
	std::cout << "Test for the Increment() method:" << std::endl;
	std::cout << counts << std::endl;

	std::cout << std::endl;
	std::cout << "Test for the GetCount() method:" << std::endl;
	std::cout << gene6 << ", " << counts.GetCount(gene6) << std::endl;

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

