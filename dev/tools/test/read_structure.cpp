/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>
#include <string>

#include "read_structure.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	char definition[] = "5C3X18U6C3X8M";
	ReadStructure structure = ReadStructure(definition);
	std::cout << structure << std::endl;

	////////////////////////////////////////////////////////////////////////////
	char read1[] = "CGCCGAAATCTTCAGCGTTCCGATAGAACCGTTCCTCTGGCA";
	std::map<char, std::string> sequences = structure.GetSequences(read1);
	std::cout << sequences['U'] << std::endl;
	std::cout << sequences['C'] << std::endl;
	std::cout << sequences['M'] << std::endl;
	std::cout << sequences['X'] << std::endl;

	////////////////////////////////////////////////////////////////////////////
	std::cout << structure.GetLength('U') << std::endl;
	std::cout << structure.GetLength('C') << std::endl;
	std::cout << structure.GetLength('M') << std::endl;
	std::cout << structure.GetLength('X') << std::endl;

	////////////////////////////////////////////////////////////////////////////
	std::cout << structure.GetMinReadLength() << std::endl;

	return 0;
}
