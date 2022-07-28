/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>
#include <seqan/sequence.h>

#include "read2.hpp"

int main(int argc, char** argv)
{
	seqan::CharString meta = "A01366:105:H5M7KDMXY:1:1101:6316:1094";
	seqan::CharString sequence = "GGCCGCACTCTCTCTGATTACAACATCCAGAAAGAGTCGACCCTGCACCTGGTCCTCCGTCTGAGGGGTGGCTATTAATTATTCGGTCTGCATTCCCAGTG";
	seqan::CharString quality = "FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFFF:FFFFFFFFFFFF::FFFFFFFFFFFFFF";

	Read2 read2 = Read2(toCString(meta), toCString(sequence), toCString(quality));
	std::cout << read2 << std::endl;

	return 0;
}
