/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>

#include "mapping.hpp"
#include "mappings.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////

	Mapping mapping1 = Mapping("GeneA", 5);
	Mapping mapping2 = Mapping("GeneA", 5);
	Mapping mapping3 = Mapping("GeneA", 6);
	Mapping mapping4 = Mapping("GeneA", 7);
	Mapping mapping5 = Mapping("GeneA", 7);
	Mapping mapping6 = Mapping("GeneB", 2);
	Mapping mapping7 = Mapping("GeneB", 9);

	Mappings mappings;
	mappings.Insert(mapping1);
	mappings.Insert(mapping2);
	mappings.Insert(mapping3);
	mappings.Insert(mapping4);
	mappings.Insert(mapping5);
	mappings.Insert(mapping6);
	mappings.Insert(mapping7);
	
	////////////////////////////////////////////////////////////////////////////
	
	std::cout << std::endl;
	std::cout << "Test to display all the mappings:" << std::endl;
	for (auto& m : mappings)
	{
		std::cout << m << ", " << m.GetScore() << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

