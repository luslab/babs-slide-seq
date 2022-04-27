/*
 * nourdine.bah@crick.ac.uk
 */

#include "utils.hpp"
#include "read_structure.hpp"
#include "read1.hpp"

#include <iostream>

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Read1::Read1()
{
}

// ----------------------------------------------------------------------------
// Full constructor
// ----------------------------------------------------------------------------

Read1::Read1(ReadStructure& structure, char* sequence, int max_distance)
{
	std::string seq(sequence);
	this->sequence = seq;

	if ( seq.size() < structure.GetMinReadLength() )
	{
		this->long_enough = false;
		this->up_primer = DummySeq( structure.GetLength('U') );
		this->up_distance = structure.GetLength('U');
		this->up_matched = false;
		this->barcode = DummySeq( structure.GetLength('C') );
		this->barcode_matched = false;
		this->umi = DummySeq( structure.GetLength('M') );
	}

	else
	{
		this->long_enough = true;
		std::map<char, std::string> sequences = structure.GetSequences(sequence);
		this->up_primer = sequences['U'];
		this->barcode = sequences['C'];
		this->umi = sequences['M'];

		this->up_distance = HammingDist(this->up_primer, "TCTTCAGCGTTCCCGAGA");
		if ( this->up_distance >= max_distance )
		{
			this->up_matched = false;
		}
		else {
			this->up_matched = true;
		}
	}
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetSeq()
// ----------------------------------------------------------------------------

std::string Read1::GetSeq()
{
	return this->sequence;
}

// ----------------------------------------------------------------------------
// GetSeqStatus()
// ----------------------------------------------------------------------------

std::string Read1::GetSeqStatus()
{
	if ( this->long_enough ) { return "LONG_ENOUGH"; }
	else { return "TOO_SHORT"; }
}

// ----------------------------------------------------------------------------
// GetUP()
// ----------------------------------------------------------------------------

std::string Read1::GetUP()
{
	return this->up_primer;
}

// ----------------------------------------------------------------------------
// GetUPDist()
// ----------------------------------------------------------------------------

int Read1::GetUPDist()
{
	return this->up_distance;
}

// ----------------------------------------------------------------------------
// GetUPStatus()
// ----------------------------------------------------------------------------

std::string Read1::GetUPStatus()
{
	if ( this->up_matched ) { return "MATCHED"; }
	else { return "UNMATCHED"; }
}

// ----------------------------------------------------------------------------
// GetBarcode()
// ----------------------------------------------------------------------------

std::string Read1::GetBarcode()
{
	return this->barcode;
}

// ----------------------------------------------------------------------------
// GetBarcodeStatus()
// ----------------------------------------------------------------------------

std::string Read1::GetBarcodeStatus()
{
	if ( this->barcode_matched ) { return "MATCHED"; }
	else { return "UNMATCHED"; }
}

// ----------------------------------------------------------------------------
// GetUMI()
// ----------------------------------------------------------------------------

std::string Read1::GetUMI()
{
	return this->umi;
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, Read1& read1)
{
	out << "(" << read1.GetBarcode() << ", " << read1.GetUMI() << ", " << read1.GetUPDist() << ")";
	return out;
}

