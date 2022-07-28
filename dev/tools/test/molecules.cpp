/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>

#include "record.hpp"
#include "molecule.hpp"
#include "molecules.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////

	Record rec1 = Record(2, "read1", 10, "GeneA");
	Record rec2 = Record(3, "read1", 14, "GeneB");
	Record rec3 = Record(4, "read2", 14, "GeneB");
	Record rec4 = Record(4, "read2", 14, "GeneB");
	Record rec5 = Record(4, "read3", 15, "GeneB");
	Record rec6 = Record(4, "read3", 14, "GeneC");

	
	Molecule mol1 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec1);
	Molecule mol2 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec2);
	Molecule mol3 = Molecule("TGCATGCATGCATGCA", "ACGTACGT", rec3);

	Molecules molecules;

	molecules.Insert(mol1);
	molecules.Insert(mol2);
	molecules.Insert(mol3);

	////////////////////////////////////////////////////////////////////////////
	
	std::cout << std::endl;
	std::cout << "Test to display all the molecules:" << std::endl;
	for (auto& mol : molecules)
	{
		std::cout << mol << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

