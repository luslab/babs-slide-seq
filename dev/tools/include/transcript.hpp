/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __TRANSCRIPT_HPP
#define __TRANSCRIPT_HPP

#include <string>

#include "read1.hpp"
#include "read2.hpp"

class Transcript
{
	private:
		Read1 read1;
		Read2 read2;
	
	public:
		Transcript(Read1, Read2);
		std::string GetMeta();
		std::string GetSeq();
		std::string GetQual();
		friend std::ostream& operator<<(std::ostream&, Transcript&);
};

#endif
