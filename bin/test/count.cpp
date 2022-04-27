/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <set>
#include <iostream>

#include "gene.hpp"
#include "count.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	
	Gene gene1 = Gene("ENSMUS000", "GeneA");
	Gene gene2 = Gene("ENSMUS001", "GeneB");
	Gene gene3 = Gene("ENSMUS002", "GeneC");
	Gene gene4 = Gene("ENSMUS002", "GeneC");
	Gene gene5 = Gene("ENSMUS002", "GeneD");

	Count count = Count(gene1);
	std::cout << std::endl;
	std::cout << "Test for the Increment() method:" << std::endl;
	std::cout << count << std::endl;
	count.Increment();
	std::cout << count << std::endl;
	count.Increment();
	count.Increment();
	std::cout << count << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::set<Count> counts;
	counts.insert( Count(gene1) );
	counts.insert( Count(gene2) );
	counts.insert( Count(gene3) );
	counts.insert( Count(gene4) );
	counts.insert( Count(gene5) );

	std::cout << std::endl;
	std::cout << "Test for Gene in set:" << std::endl;
	for (auto& count : counts)
	{
		std::cout << count << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

