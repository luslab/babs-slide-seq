/*
 * nourdine.bah@crick.ac.uk
 */

#include <tuple>
#include <set>
#include <iostream>

#include "mapping.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	
	Mapping mapping1 = Mapping("GeneA", 5);
	Mapping mapping2 = Mapping("GeneA", 5);
	Mapping mapping3 = Mapping("GeneA", 6);
	Mapping mapping4 = Mapping("GeneA", 7);
	Mapping mapping5 = Mapping("GeneB", 2);
	Mapping mapping6 = Mapping("GeneB", 9);
	
	std::cout << std::endl;
	std::cout << "Test of the ostream" << std::endl;
	std::cout << mapping1 << std::endl;
	std::cout << mapping2 << std::endl;
	std::cout << mapping3 << std::endl;
	std::cout << mapping4 << std::endl;
	std::cout << mapping5 << std::endl;
	std::cout << mapping6 << std::endl;

	////////////////////////////////////////////////////////////////////////////
	
	std::cout << std::endl;
	std::cout << "Test of the copy constructor:" << std::endl;
	Mapping copied_mapping = Mapping(mapping1);
	std::cout << mapping1 << std::endl;
	std::cout << copied_mapping << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the + operator with same genes:" << std::endl;
	Mapping m_mapping1;
	try {
		m_mapping1 = mapping1 + mapping2 + mapping3;
	}
	catch(std::tuple<Mapping, Mapping> mappings)
	{
		std::cerr << "It should not fail here" << std::endl;
	}
	std::cout << mapping1 << ", " << mapping2 << ", " << mapping3 << ", " << m_mapping1 << std::endl;

	std::cout << std::endl;
	std::cout << "Test of the + operator with different genes:" << std::endl;
	Mapping m_mapping2;
	try {
		m_mapping1 = mapping1 + mapping5;
	}
	catch(std::tuple<Mapping, Mapping> mappings)
	{
		std::cerr << "It should fail here" << std::endl;
	}
	std::cout << mapping1 << ", " << mapping5 << ", " << m_mapping2 << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the GetScore() method:" << std::endl;
	std::cout << m_mapping1 << ", " << m_mapping1.GetScore() << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test putting mappings in a set:" << std::endl;
	std::set<Mapping> mappings;
	mappings.insert(mapping1);
	mappings.insert(mapping2);
	mappings.insert(mapping3);
	mappings.insert(mapping4);
	mappings.insert(mapping5);
	mappings.insert(mapping6);
	for (auto& m : mappings)
	{
		std::cout << m << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

