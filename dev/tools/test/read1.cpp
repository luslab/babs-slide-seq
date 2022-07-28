/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>
#include <string>

#include "read1.hpp"

int main(int argc, char** argv)
{
	char definition[] = "5C3X18U6C3X8M";
	ReadStructure structure = ReadStructure(definition);
	std::cout << structure << std::endl;

	char sequence[] = "CGCCGAAATCTTCAGCGTTCCGATAGAACCGTTCCTCTGGCA";
	std::cout << structure.GetStruct() << std::endl;
	std::cout << sequence << std::endl;

	Read1 read1 = Read1(structure, sequence, 3);
	std::cout << read1 << std::endl;

	return 0;
}
