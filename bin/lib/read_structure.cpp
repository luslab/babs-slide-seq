/*
 * nourdine.bah@crick.ac.uk
 */

#include <vector>
#include <map>
#include <fstream>

#include "read_structure.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

ReadStructure::ReadStructure()
{
	this->definition = std::string("8C18U6C3X8M");
	this->structure = ReadStructure::StructFromDef(this->definition);
	this->lengths = this->GetLengths();
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

ReadStructure::ReadStructure(char* definition)
{
	this->definition = std::string(definition);
	this->structure = ReadStructure::StructFromDef(this->definition);
	this->lengths = this->GetLengths();
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetDef()
// ----------------------------------------------------------------------------

std::string ReadStructure::GetDef()
{
	return definition;
}

// ----------------------------------------------------------------------------
// GetStruct()
// ----------------------------------------------------------------------------

std::string ReadStructure::GetStruct()
{
	return structure;
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// StructFromDef()
// ----------------------------------------------------------------------------

std::string ReadStructure::StructFromDef(std::string definition)
{
	int size = 0;
	std::string size_str = "";
	std::vector< std::tuple<char, int> > segments;

	// parsing the numbers
	for (char const &c: definition)
	{
		if ( isdigit(c) )
		{
			size_str += c;
		}
		else
		{
			size = atoi( size_str.c_str() );
			size_str = "";
			segments.push_back( std::tuple<int, char>{size, c} );
		}
	}

	// create structure
	std::string structure = "";
	for (auto& segment: segments)
	{
		for (int i=0; i<std::get<0>(segment); i++)
		{
			structure.push_back(std::get<1>(segment));
		}
	}

	return structure;
}

// ----------------------------------------------------------------------------
// GetSequences()
// ----------------------------------------------------------------------------

std::map<char, std::string> ReadStructure::GetSequences(char* sequence)
{
	std::map<char, std::string> sequences;

	for (int i=0; i<this->structure.size(); i++)
	{
		sequences[ this->structure[i] ].push_back( sequence[i] );
	}

	return sequences;
}

// ----------------------------------------------------------------------------
// GetLengths()
// ----------------------------------------------------------------------------

std::map<char, int> ReadStructure::GetLengths()
{
	std::map<char, int> lengths;

	for (int i=0; i<this->structure.size(); i++)
	{
		lengths[ this->structure[i] ]++;
	}

	return lengths;
}

// ----------------------------------------------------------------------------
// GetLength()
// ----------------------------------------------------------------------------

int ReadStructure::GetLength(char c)
{
	return this->lengths[c];
}

// ----------------------------------------------------------------------------
// GetMinReadLength()
// ----------------------------------------------------------------------------

int ReadStructure::GetMinReadLength()
{
	return this->structure.size();
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, ReadStructure& read_structure)
{
	out << "(" << read_structure.GetDef() << ", " << read_structure.GetStruct() << ")";
	return out;
}

