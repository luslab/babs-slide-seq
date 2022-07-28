/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __COUNTER_HPP
#define __COUNTER_HPP

#include <string>
#include <set>
#include <vector>
#include <map>

#include "gene.hpp"
#include "counts.hpp"
#include "counter.hpp"

class Counter
{
	private:
		std::set<Gene> genes;
		std::vector<Gene> genes_v;
		std::map<std::string, std::string> symbols;
		std::map<Gene, unsigned long long> positions;
		std::map<std::string, Counts> counts;
	
	public:
		Counter();
		Counter(std::set<Gene>);
		std::vector<Gene> GetGenes() const;
		std::vector<std::string> GetBarcodes() const;
		void Increment(char*, char*);
		unsigned long long GetNumberOfGenes() const;
		unsigned long long GetNumberOfBarcodes() const;
		unsigned long long GetNumberOfEntries() const;
		void GetDGE(std::ostream&, std::ostream&, std::ostream&);
		friend std::ostream& operator<<(std::ostream&, const Counter&);
};

#endif

