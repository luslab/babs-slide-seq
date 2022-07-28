/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <tuple>
#include <vector>
#include <fstream>
#include <algorithm>

#include "mapping.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Mapping::Mapping()
{
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Mapping::Mapping(char* gene, int score)
{
	this->gene = std::string(gene);
	this->scores = {score};
}

// ----------------------------------------------------------------------------
// Alternative constructor
// ----------------------------------------------------------------------------

Mapping::Mapping(std::string gene, int score)
{
	this->gene = gene;
	this->scores = {score};
}

// ----------------------------------------------------------------------------
// Copy constructor
// ----------------------------------------------------------------------------

Mapping::Mapping(const Mapping& mapping)
{
	this->gene = mapping.GetGene();
	this->scores= std::vector<int>();
	for (auto& score : mapping.scores) // the access modifiers work on class level, and not on object level!
	{
		this->scores.push_back(score);
	}
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetGene()
// ----------------------------------------------------------------------------

std::string Mapping::GetGene() const
{
	return gene;
}

// ----------------------------------------------------------------------------
// GetScores()
// ----------------------------------------------------------------------------

std::vector<int> Mapping::GetScores() const
{
	return scores;
}

// ----------------------------------------------------------------------------
// GetScore()
// ----------------------------------------------------------------------------

int Mapping::GetScore() const
{
	std::vector<int>::const_iterator it;
	it = std::max_element(scores.begin(), scores.end());

	return *it;
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// GetScoreString()
// ----------------------------------------------------------------------------

std::string Mapping::GetScoreString() const
{
	unsigned long long i = 0;
	std::string str = "";

	for (auto& score : scores)
	{
		str.append( std::to_string(score) );

		i++;

		if ( i < scores.size() )
		{
			str.append(":");
		}
	}

	return str;
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Mapping& Mapping)
{
	out << "(" << Mapping.GetGene() << ", " << Mapping.GetScoreString() << ")";
	return out;
}

// ----------------------------------------------------------------------------
// plus
// ----------------------------------------------------------------------------

Mapping operator+(const Mapping& mapping1, const Mapping& mapping2)
{
	if ( mapping1.GetGene() ==  mapping2.GetGene() )
	{
		Mapping mapping = Mapping(mapping1);

		for (auto& score : mapping2)
		{
			mapping.scores.push_back(score);
		}

		return mapping;
	}
	else
	{
		throw( std::make_tuple(mapping1, mapping2) );
	}
}

// ----------------------------------------------------------------------------
// lesser than
// ----------------------------------------------------------------------------

bool operator<(const Mapping& Mapping1, const Mapping& Mapping2)
{
	return Mapping1.GetGene() < Mapping2.GetGene();
}

// ============================================================================
// Iterators
// ============================================================================

// ----------------------------------------------------------------------------
// begin()
// ----------------------------------------------------------------------------

std::vector<int>::const_iterator Mapping::begin() const
{
	return scores.begin();
}

// ----------------------------------------------------------------------------
// end()
// ----------------------------------------------------------------------------

std::vector<int>::const_iterator Mapping::end() const
{
	return scores.end();
}

