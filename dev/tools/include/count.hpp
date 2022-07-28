/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __COUNT_HPP
#define __COUNT_HPP

#include "gene.hpp"

class Count
{
	private:
		Gene gene;
		unsigned long long count;
	
	public:
		Count();
		Count(Gene);
		Gene GetGene() const;
		unsigned long long GetCount() const;
		void Increment();
		friend std::ostream& operator<<(std::ostream&, const Count&);
		friend bool operator<(const Count&, const Count&);
};

#endif

