/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <tuple>
#include <set>
#include <iostream>

#include "mapping.hpp"
#include "mappings.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Mappings::Mappings()
{
	this->mappings = std::set<Mapping>();
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Mappings::Mappings(Mapping mapping)
{
	this->mappings = std::set<Mapping>();
	this->mappings.insert(mapping);
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetSize()
// ----------------------------------------------------------------------------

size_t Mappings::GetSize() const
{
	return mappings.size();
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// Insert()
// ----------------------------------------------------------------------------

void Mappings::Insert(Mapping mapping)
{
	std::set<Mapping>::iterator it;
	it = mappings.find(mapping);
	if ( it == mappings.end() )
	{
		mappings.insert(mapping);
	}
	else
	{
		Mapping new_map;

		try
		{
			new_map = *it + mapping;
		}
		catch(std::tuple<Mapping, Mapping> maps)
		{
			std::cerr << "Trying to merge two mappings with different genes" << std::endl;
			std::cerr << std::get<0>(maps) << std::endl;
			std::cerr << std::get<1>(maps) << std::endl;
			exit(EXIT_FAILURE);
		}

		mappings.erase(it);
		mappings.insert(new_map);
	}
}

// ============================================================================
// Iterators
// ============================================================================

// ----------------------------------------------------------------------------
// begin()
// ----------------------------------------------------------------------------

std::set<Mapping>::iterator Mappings::begin() const
{
	return mappings.begin();
}

// ----------------------------------------------------------------------------
// end()
// ----------------------------------------------------------------------------

std::set<Mapping>::iterator Mappings::end() const
{
	return mappings.end();
}

