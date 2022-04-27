/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __COUNT_HPP
#define __COUNT_HPP

#include <map>

#include "gene.hpp"

class Counts
{
	private:
		std::map<Gene, unsigned long long> counts;
	
	public:
		Counts();
		Counts(Gene);
		unsigned long long GetSize() const;
		void Increment(Gene);
		unsigned long long GetCount(Gene);
		std::map<Gene, unsigned long long> GetCounts();
		friend std::ostream& operator<<(std::ostream&, const Counts&);
		std::map<Gene, unsigned long long>::const_iterator begin() const;
		std::map<Gene, unsigned long long>::const_iterator end() const;
};

#endif

