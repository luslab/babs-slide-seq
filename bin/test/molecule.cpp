/*
 * nourdine.bah@crick.ac.uk
 */

#include <tuple>
#include <iostream>

#include "molecule.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////

	Record rec1 = Record(2, "read1", 10, "GeneA");
	Record rec2 = Record(3, "read1", 14, "GeneB");
	Record rec3 = Record(4, "read2", 14, "GeneB");
	Record rec4 = Record(4, "read2", 14, "GeneB");
	Record rec5 = Record(4, "read3", 15, "GeneB");
	Record rec6 = Record(4, "read3", 14, "GeneC");
	Record rec7 = Record(5, "read3", 14, "GeneC");
	Record rec8 = Record(6, "read3", 18, "GeneC");
	Record rec9 = Record(7, "read3", 20, "GeneD");

	////////////////////////////////////////////////////////////////////////////

	Molecule molecule = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec1);
	molecule.Insert(rec2);
	molecule.Insert(rec3);
	molecule.Insert(rec4);
	molecule.Insert(rec5);
	molecule.Insert(rec6);
	molecule.Insert(rec7);
	molecule.Insert(rec8);
	molecule.Insert(rec9);

	std::cout << std::endl;
	std::cout << "Test of the ostream:" << std::endl;
	std::cout << molecule << std::endl;

	std::cout << std::endl;
	std::cout << "Test to display all the records:" << std::endl;
	for (auto& rec : molecule)
	{
		std::cout << rec << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Test for the GetGenes() method:" << std::endl;
	for (auto& gene : molecule.GetGenes())
	{
		std::cout << gene << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Test for the GetGeneString() method:" << std::endl;
	std::cout << molecule.GetGeneString() << std::endl;

	////////////////////////////////////////////////////////////////////////////
	
	std::cout << std::endl;
	std::cout << "Test of the copy constructor:" << std::endl;
	Molecule mol = Molecule(molecule);
	std::cout << mol << std::endl;
	std::cout << molecule << std::endl;

	////////////////////////////////////////////////////////////////////////////
	
	Molecule mol1 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec1);
	Molecule mol2 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec2);
	Molecule mol3 = Molecule("TGCATGCATGCATGCA", "ACGTACGT", rec3);

	std::cout << std::endl;
	std::cout << "Test of the + operator with same barcodes:" << std::endl;
	Molecule m_mol1;
	try {
		m_mol1 = mol1 + mol2;
	}
	catch(std::tuple<Molecule, Molecule> molecules)
	{
		std::cerr << "It should not fail here" << std::endl;
	}
	std::cout << mol1 << ", " << mol2 << ", " << m_mol1 << std::endl;

	std::cout << std::endl;
	std::cout << "Test of the + operator with different barcodes:" << std::endl;
	Molecule m_mol2;
	try
	{
		m_mol2 = mol2 + mol3;
	}
	catch(std::tuple<Molecule, Molecule> molecules)
	{
		std::cerr << "It should fail here" << std::endl;
		std::cerr << "Trying to merge two molecules with different (barcode, UMI) pairs:" << std::endl;
		std::cerr << std::get<0>(molecules) << std::endl;
		std::cerr << std::get<1>(molecules) << std::endl;
	}
	std::cout << mol2 << ", " << mol3 << ", " << m_mol2 << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the ExtractMappings() and GetMappings() methods:" << std::endl;
	Mappings mappings = molecule.GetMappings(true);
	for (auto& m : mappings)
	{
		std::cout << m << ", " << m.GetScore() << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Test of the GetMaxScore() method:" << std::endl;
	std::cout << molecule.GetMaxScore() << std::endl;

	////////////////////////////////////////////////////////////////////////////

	Record r1 = Record(1, "read1", 1, "GeneA");
	Record r2 = Record(2, "read2", 2, "GeneB");
	Record r3 = Record(3, "read3", 3, "GeneC");
	Record r4 = Record(4, "read4", 1, "GeneD");
	Record r5 = Record(5, "read5", 2, "GeneE");
	Record r6 = Record(6, "read6", 3, "GeneF");
	Record r7 = Record(7, "read7", 1, "GeneG");
	Record r8 = Record(8, "read8", 2, "GeneH");
	Record r9 = Record(9, "read9", 3, "GeneI");

	Molecule mm1 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", r1);
	mm1.Insert(r2);
	mm1.Insert(r3);
	mm1.Insert(r4);
	mm1.Insert(r5);

	Molecule mm2 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", r3);
	mm2.Insert(r6);
	mm2.Insert(r8);

	mm1.ExtractMappings();
	mm2.ExtractMappings();

	std::cout << std::endl;
	std::cout << "Test of the IsThereAMaxima() method:" << std::endl;
	std::cout << mm1 << ", this should be true: " << mm1.IsThereAMaxima() << std::endl;
	std::cout << mm2 << ", this should be false :" << mm2.IsThereAMaxima() << std::endl;

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the GetRecordTags() method:" << std::endl;

	for (auto& rec : mm1)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm1.GetRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	for (auto& rec : mm2)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm2.GetRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Molecule mm3 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", r1);
	mm3.ExtractMappings();
	for (auto& rec : mm3)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm3.GetRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Record rr1 = Record(1, "read1", 1, "GeneA");
	Record rr2 = Record(2, "read2", 2, "GeneA");
	Record rr3 = Record(3, "read2", 2, "GeneA");
	Molecule mm4 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rr1);
	mm4.Insert(rr2);
	mm4.Insert(rr3);
	mm4.ExtractMappings();
	for (auto& rec : mm4)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm4.GetRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Molecule mm5 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rr2);
	mm5.Insert(rr3);
	mm5.ExtractMappings();
	for (auto& rec : mm5)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag] : mm5.GetRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	std::cout << std::endl;
	std::cout << "Test of the GetFrequencyBasedRecordTags() method:" << std::endl;

	mm1.ComputeFrequencies();
	for (auto& rec : mm1)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm1.GetFrequencyBasedRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	mm2.ComputeFrequencies();
	for (auto& rec : mm2)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag]: mm2.GetFrequencyBasedRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Molecule mm6 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec1);
	mm6.ComputeFrequencies();
	for (auto& rec : mm6)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag] : mm6.GetFrequencyBasedRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Molecule mm7 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec1);
	mm7.Insert(rec3);
	mm7.Insert(rec4);
	mm7.Insert(rec5);
	mm7.Insert(rec6);
	mm7.Insert(rec7);
	mm7.Insert(rec8);
	mm7.Insert(rec9);
	mm7.ComputeFrequencies();
	for (auto& rec : mm7)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag] : mm7.GetFrequencyBasedRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	Molecule mm8 = Molecule("ACGTACGTACGTACGT", "ACGTACGT", rec2);
	mm8.Insert(rec3);
	mm8.Insert(rec6);
	mm8.Insert(rec7);
	mm8.ComputeFrequencies();
	for (auto& rec : mm8)
	{
		std::cout << rec << ", " << rec.GetScore() << std::endl;
	}
	for (auto& [pos, tag] : mm8.GetFrequencyBasedRecordTags())
	{
		std::cout << pos << ", " << tag << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////

	return 0;
}

