/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __GENE_HPP
#define __GENE_HPP

#include <string>

class Gene
{
	private:
		std::string id;
		std::string symbol;
	
	public:
		Gene();
		Gene(char*, char*);
		Gene(std::string, std::string);
		std::string GetID() const;
		std::string GetSymbol() const;
		friend std::ostream& operator<<(std::ostream&, const Gene&);
		friend bool operator<(const Gene&, const Gene&);
};

#endif

