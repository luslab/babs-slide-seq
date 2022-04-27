/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __READ2_HPP
#define __READ2_HPP

#include <seqan/basic.h>
#include <seqan/sequence.h>

class Read2
{
	private:
		std::string meta;
		std::string sequence;
		std::string quality;
	
	public:
		Read2();
		Read2(char*, char*, char*);
		std::string GetMeta();
		std::string GetSeq();
		std::string GetQual();
		friend std::ostream& operator<<(std::ostream&, Read2&);
};

#endif
