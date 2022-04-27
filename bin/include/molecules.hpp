/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __MOLECULES_HPP
#define __MOLECULES_HPP

#include <set>

#include "molecule.hpp"

class Molecules
{
	private:
		std::set<Molecule> molecules;
	
	public:
		Molecules();
		Molecules(Molecule);
		size_t GetSize() const;
		void Insert(Molecule);
		std::set<Molecule>::iterator begin();
		std::set<Molecule>::iterator end();
};

#endif

