/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __READ1_HPP
#define __READ1_HPP

#include "read_structure.hpp"

class Read1
{
	private:
		std::string sequence;
		bool long_enough;

		std::string up_primer;
		int up_distance;
		bool up_matched;

		std::string barcode;
		bool barcode_matched;

		std::string umi;
	
	public:
		Read1();
		Read1(ReadStructure&, char*, int);
		std::string GetSeq();
		std::string GetSeqStatus();
		std::string GetUP();
		std::string GetUPStatus();
		int GetUPDist();
		std::string GetBarcode();
		std::string GetBarcodeStatus();
		std::string GetUMI();
		friend std::ostream& operator<<(std::ostream&, Read1&);
};

#endif
