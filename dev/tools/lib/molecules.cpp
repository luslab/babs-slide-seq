/*
 * nourdine.bah@crick.ac.uk
 */

#include <tuple>
#include <set>
#include <iostream>
#include <fstream>

#include "molecule.hpp"
#include "molecules.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Molecules::Molecules()
{
	this->molecules = std::set<Molecule>();
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Molecules::Molecules(Molecule molecule)
{
	this->molecules = std::set<Molecule>();
	this->molecules.insert(molecule);
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetSize()
// ----------------------------------------------------------------------------

size_t Molecules::GetSize() const
{
	return molecules.size();
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// Insert()
// ----------------------------------------------------------------------------

void Molecules::Insert(Molecule molecule)
{
	std::set<Molecule>::iterator it;
	it = molecules.find(molecule);
	if ( it == molecules.end() )
	{
		molecules.insert(molecule);
	}
	else
	{
		Molecule new_mol(*it);

		try
		{
			new_mol = new_mol + molecule;
		}
		catch(std::tuple<Molecule, Molecule> mols)
		{
			std::cerr << "Trying to merge two molecules with different (barcode, UMI) pairs" << std::endl;
			std::cerr << std::get<0>(mols) << std::endl;
			std::cerr << std::get<1>(mols) << std::endl;
			exit(EXIT_FAILURE);
		}

		molecules.erase(it);
		molecules.insert(new_mol);
	}
}

// ============================================================================
// Iterators
// ============================================================================

// ----------------------------------------------------------------------------
// begin()
// ----------------------------------------------------------------------------

std::set<Molecule>::iterator Molecules::begin()
{
	return molecules.begin();
}

// ----------------------------------------------------------------------------
// end()
// ----------------------------------------------------------------------------

std::set<Molecule>::iterator Molecules::end()
{
	return molecules.end();
}

