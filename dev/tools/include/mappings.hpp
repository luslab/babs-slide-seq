/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __MAPPINGS_HPP
#define __MAPPINGS_HPP

#include <set>

#include "mapping.hpp"

class Mappings
{
	private:
		std::set<Mapping> mappings;
	
	public:
		Mappings();
		Mappings(Mapping);
		size_t GetSize() const;
		void Insert(Mapping);
		std::set<Mapping>::iterator begin() const;
		std::set<Mapping>::iterator end() const;
};

#endif

