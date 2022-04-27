/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __READ_STRUCTURE_HPP
#define __READ_STRUCTURE_HPP

#include <map>
#include <string>

class ReadStructure
{
	private:
		std::string	definition;
		std::string	structure;
		std::map<char, int> lengths;
	
	public:
		ReadStructure();
		ReadStructure(char*);
		std::string StructFromDef(std::string);
		std::string GetDef();
		std::string GetStruct();
		std::map<char, std::string> GetSequences(char*);
		std::map<char, int> GetLengths();
		int GetLength(char);
		int GetMinReadLength();
		friend std::ostream& operator<<(std::ostream&, ReadStructure&);
};

#endif
