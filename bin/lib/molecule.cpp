/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <tuple>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

#include "record.hpp"
#include "mapping.hpp"
#include "mappings.hpp"
#include "molecule.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Molecule::Molecule()
{
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Molecule::Molecule(char* barcode, char* umi, Record record)
{
	this->barcode = std::string(barcode);
	this->umi = std::string(umi);
	this->sequence = this->barcode + this->umi;
	this->records = {record};
}

// ----------------------------------------------------------------------------
// Copy constructor
// ----------------------------------------------------------------------------

Molecule::Molecule(const Molecule& molecule)
{
	this->barcode = molecule.GetBarcode();
	this->umi = molecule.GetUMI();
	this->sequence = molecule.GetSeq();
	this->records = std::set<Record>();
	for (auto& rec : molecule.records) // the access modifiers work on class level, and not on object level!
	{
		this->records.insert(rec);
	}
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetBarcode()
// ----------------------------------------------------------------------------

std::string Molecule::GetBarcode() const
{
	return barcode;
}

// ----------------------------------------------------------------------------
// GetUMI()
// ----------------------------------------------------------------------------

std::string Molecule::GetUMI() const
{
	return umi;
}

// ----------------------------------------------------------------------------
// GetSeq()
// ----------------------------------------------------------------------------

std::string Molecule::GetSeq() const
{
	return sequence;
}

// ----------------------------------------------------------------------------
// GetRecords()
// ----------------------------------------------------------------------------

std::set<Record> Molecule::GetRecords() const
{
	return records;
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// GetRecordString()
// ----------------------------------------------------------------------------

std::string Molecule::GetRecordString() const
{
	unsigned long long i = 0;
	std::string str = "";

	for (auto& rec: records)
	{
		str.append( std::to_string(rec.GetPos()) );

		i++;

		if ( i < records.size() )
		{
			str.append(":");
		}
	}

	return str;
}

// ----------------------------------------------------------------------------
// Insert()
// ----------------------------------------------------------------------------

void Molecule::Insert(Record record)
{
	records.insert(record);
}

// ----------------------------------------------------------------------------
// GetGenes()
// ----------------------------------------------------------------------------

std::set<std::string> Molecule::GetGenes() const
{
	std::set<std::string> genes;
	for (auto& rec : records)
	{
		genes.insert( rec.GetGene() );
	}
	return genes;
}

// ----------------------------------------------------------------------------
// GetGeneString()
// ----------------------------------------------------------------------------

std::string Molecule::GetGeneString() const
{
	unsigned long long i = 0;
	std::string str = "";
	std::set<std::string> genes = GetGenes();

	for (auto& gene: genes)
	{
		str.append(gene);

		i++;

		if ( i < genes.size() )
		{
			str.append(":");
		}
	}

	return str;
}

// ----------------------------------------------------------------------------
// ExtractMappings()
// ----------------------------------------------------------------------------

void Molecule::ExtractMappings()
{
	Mappings mappings;

	for (auto& rec : records)
	{
		Mapping mapping = Mapping(rec.GetGene(), rec.GetScore());
		mappings.Insert(mapping);
	}

	this->mappings = mappings;
}

// ----------------------------------------------------------------------------
// GetMappings()
// ----------------------------------------------------------------------------

Mappings Molecule::GetMappings(bool extract=false)
{
	if (extract)
	{
		ExtractMappings();
	}

	return mappings;
}

// ----------------------------------------------------------------------------
// GetMaxScore()
// ----------------------------------------------------------------------------

long Molecule::GetMaxScore() const
{
	long max = -1;

	for (auto& mapping : mappings)
	{
		if ( mapping.GetScore() > max )
		{
			max = mapping.GetScore();
		}
	}

	return max;
}

// ----------------------------------------------------------------------------
// IsThereAMaxima()
// ----------------------------------------------------------------------------

bool Molecule::IsThereAMaxima() const
{
	int max = GetMaxScore();
	long n = 0;

	for (auto& mapping : mappings)
	{
		if ( mapping.GetScore() == max )
		{
			n++;
		}
	}

	if ( n == 0 )
	{
		std::cerr << "Error: the maximal score value can't be found in the " << *this << " molecule" << std::endl;
		exit(EXIT_FAILURE);
	}
	else if ( n == 1 )
	{
		return true;
	}
	else
	{
		return false;
	}
}

// ----------------------------------------------------------------------------
// GetRecordTags()
// ----------------------------------------------------------------------------

std::map<unsigned long long, std::string> Molecule::GetRecordTags() const
{
	std::map<unsigned long long, std::string> map;

	if ( GetGenes().size() == 1 && records.size() == 1 )
	{
		for (auto& rec : records)
		{
			map[rec.GetPos()] = "UNIQUE";
		}

		return map;
	}

	if ( IsThereAMaxima() )
	{
		bool found = false;
		long max = GetMaxScore();

		for (auto& rec : records)
		{
			if ( !found && rec.GetScore() == max )
			{
				map[rec.GetPos()] = "INCLUDED";
				found = true;
			}
			else
			{
				map[rec.GetPos()] = "EXCLUDED";
			}
		}
	}

	else
	{
		for (auto& rec : records)
		{
			map[rec.GetPos()] = "UNRESOLVED";
		}
	}

	return map;
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Molecule& molecule)
{
	out << "(" << molecule.GetBarcode() << ", " << molecule.GetUMI() << ", " << molecule.GetRecordString() << ")";
	return out;
}

// ----------------------------------------------------------------------------
// lesser than
// ----------------------------------------------------------------------------

bool operator<(const Molecule& molecule1, const Molecule& molecule2)
{
	return molecule1.GetSeq() < molecule2.GetSeq();
}

// ----------------------------------------------------------------------------
// plus
// ----------------------------------------------------------------------------

Molecule operator+(Molecule& molecule1, Molecule& molecule2)
{
	if ( molecule1.GetSeq() ==  molecule2.GetSeq() )
	{
		Molecule molecule = Molecule(molecule1);

		for (auto& rec : molecule2)
		{
			molecule.Insert(rec);
		}

		return molecule;
	}
	else
	{
		throw( std::make_tuple(molecule1, molecule2) );
	}
}

// ============================================================================
// Iterators
// ============================================================================

// ----------------------------------------------------------------------------
// begin()
// ----------------------------------------------------------------------------

std::set<Record>::iterator Molecule::begin()
{
	return records.begin();
}

// ----------------------------------------------------------------------------
// end()
// ----------------------------------------------------------------------------

std::set<Record>::iterator Molecule::end()
{
	return records.end();
}

