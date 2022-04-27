/*
 * nourdine.bah@crick.ac.uk
 */

#include "transcript.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Transcript::Transcript(Read1 read1, Read2 read2)
{
	this->read1 = read1;
	this->read2 = read2;
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetMeta()
// ----------------------------------------------------------------------------

std::string Transcript::GetMeta()
{
	std::string meta;

	meta.append(read1.GetSeqStatus());
	meta.push_back(':');
	meta.append(read1.GetUP());
	meta.push_back(':');
	meta.append(std::to_string(read1.GetUPDist()));
	meta.push_back(':');
	meta.append(read1.GetUPStatus());
	meta.push_back(':');
	meta.append(read1.GetBarcode());
	meta.push_back(':');
	meta.append(read1.GetUMI());
	meta.push_back(':');
	meta.append(read2.GetMeta());

	return meta;
}

// ----------------------------------------------------------------------------
// GetSeq()
// ----------------------------------------------------------------------------

std::string Transcript::GetSeq()
{
	std::string sequence = this->read2.GetSeq();
	return sequence;
}

// ----------------------------------------------------------------------------
// GetQual()
// ----------------------------------------------------------------------------

std::string Transcript::GetQual()
{
	std::string quality = this->read2.GetQual();
	return quality;
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, Transcript& transcript)
{
	out << "(" << transcript.GetMeta() << ")";
	return out;
}

