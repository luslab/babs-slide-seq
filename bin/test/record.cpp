/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>
#include <set>

#include "record.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	Record rec1 = Record(2, "read1", 10, "GeneA");
	Record rec2 = Record(3, "read1", 14, "GeneB");
	Record rec3 = Record(4, "read2", 14, "GeneB");
	std::cout << rec1 << std::endl;
	std::cout << rec2 << std::endl;
	std::cout << rec3 << std::endl;

	////////////////////////////////////////////////////////////////////////////
	std::set<Record> records;
	records.insert(rec1);
	records.insert(rec1);
	records.insert(rec2);
	records.insert(rec2);
	records.insert(rec3);
	records.insert(rec3);
	for (auto& rec: records)
	{
		std::cout << "In set, " << rec << std::endl;
	}

	return 0;
}
