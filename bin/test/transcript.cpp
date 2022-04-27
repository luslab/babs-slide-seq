/*
 * nourdine.bah@crick.ac.uk
 */

#include "transcript.hpp"

int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	char definition[] = "5C3X18U6C3X8M";
	ReadStructure structure = ReadStructure(definition);
	std::cout << structure << std::endl;

	////////////////////////////////////////////////////////////////////////////
	char sequence1[] = "CGCCGAAATCTTCAGCGTTCCGATAGAACCGTTCCTCTGGCA";
	std::cout << structure.GetStruct() << std::endl;
	std::cout << sequence1 << std::endl;

	////////////////////////////////////////////////////////////////////////////
	Read1 read1 = Read1(structure, sequence1, 3);
	std::cout << read1 << std::endl;

	////////////////////////////////////////////////////////////////////////////
	seqan::CharString meta = "A01366:105:H5M7KDMXY:1:1101:6316:1094";
	seqan::CharString sequence2 = "GGCCGCACTCTCTCTGATTACAACATCCAGAAAGAGTCGACCCTGCACCTGGTCCTCCGTCTGAGGGGTGGCTATTAATTATTCGGTCTGCATTCCCAGTG";
	seqan::CharString quality = "FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFF::FFFFFFFFFFFFFF";

	////////////////////////////////////////////////////////////////////////////
	Read2 read2 = Read2(toCString(meta), toCString(sequence2), toCString(quality));
	std::cout << read2 << std::endl;

	////////////////////////////////////////////////////////////////////////////
	Transcript transcript = Transcript(read1, read2);
	std::cout << transcript << std::endl;
	std::cout << transcript.GetMeta() << std::endl;
	std::cout << transcript.GetSeq() << std::endl;
	std::cout << transcript.GetQual() << std::endl;

	return 0;
}
