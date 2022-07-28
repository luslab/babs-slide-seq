/*
 * nourdine.bah@crick.ac.uk
 */

#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "gene.hpp"
#include "counts.hpp"
#include "counter.hpp"

// ============================================================================
// Constructors
// ============================================================================

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------

Counter::Counter()
{
}

// ----------------------------------------------------------------------------
// Minimal constructor
// ----------------------------------------------------------------------------

Counter::Counter(std::set<Gene> genes)
{
	// The genes
	this->genes = genes;
	this->genes_v = std::vector<Gene>(this->genes.size());
	std::copy(this->genes.begin(), this->genes.end(), this->genes_v.begin());
	std::sort(this->genes_v.begin(), this->genes_v.end(), [](Gene g1, Gene g2) { return g1.GetSymbol() < g2.GetSymbol(); });

	// Maps gene IDs to gene symbols
	this->symbols = std::map<std::string, std::string>();
	for (auto& gene : genes)
	{
		this->symbols[ gene.GetID() ] = gene.GetSymbol();
	}

	// Maps ordered genes to their positions
	this->positions = std::map<Gene, unsigned long long>();
	for (unsigned long long i=0; i<this->genes_v.size(); i++)
	{
		this->positions[ genes_v[i] ] = i + 1;
	}

	// Counts
	this->counts = std::map<std::string, Counts>();
}

// ============================================================================
// Getters
// ============================================================================

// ----------------------------------------------------------------------------
// GetGenes()
// ----------------------------------------------------------------------------

std::vector<Gene> Counter::GetGenes() const
{
	return genes_v;
}

// ----------------------------------------------------------------------------
// GetBarcodes()
// ----------------------------------------------------------------------------

std::vector<std::string> Counter::GetBarcodes() const
{
	std::vector<std::string> barcodes;

	for (auto& [bcd, cnts] : counts)
	{
		barcodes.push_back(bcd);
	}

	std::sort(barcodes.begin(), barcodes.end(), [](std::string bcd1, std::string bcd2) { return bcd1 < bcd2; });

	return barcodes;
}

// ============================================================================
// Methods
// ============================================================================

// ----------------------------------------------------------------------------
// Increment()
// ----------------------------------------------------------------------------

void Counter::Increment(char* barcode, char* id)
{
	std::string barcode_str(barcode);
	std::string id_str(id);
	Gene gene = Gene(id, symbols[id]);
	counts[barcode_str].Increment(gene);
}

// ----------------------------------------------------------------------------
// GetNumberOfGenes()
// ----------------------------------------------------------------------------

unsigned long long Counter::GetNumberOfGenes() const
{
	return (unsigned long long)genes.size();
}

// ----------------------------------------------------------------------------
// GetNumberOfBarcodes()
// ----------------------------------------------------------------------------

unsigned long long Counter::GetNumberOfBarcodes() const
{
	return (unsigned long long)counts.size();
}

// ----------------------------------------------------------------------------
// GetNumberOfEntries()
// ----------------------------------------------------------------------------

unsigned long long Counter::GetNumberOfEntries() const
{
	unsigned long long n = 0;
	for (auto& [bcd, cnts] : counts)
	{
		n += cnts.GetSize();
	}
	return n;
}

// ----------------------------------------------------------------------------
// GetDGE()
// ----------------------------------------------------------------------------

void Counter::GetDGE(std::ostream& feat, std::ostream& bcd, std::ostream& mtx)
{
	// Barcodes
	std::vector<std::string> barcodes = GetBarcodes();

	// Barcodes positions
	std::map<std::string, unsigned long long> pos;
	for (unsigned long long i=0; i<barcodes.size(); i++)
	{
		pos[ barcodes[i] ] = i + 1;
	}

	// Barcode sequences
	for (auto& barcode : barcodes)
	{
		bcd << barcode << std::endl;
	}

	// Genes
	for (auto& gene : genes_v)
	{
		feat << gene.GetID() << "\t" << gene.GetSymbol() << "\tGene Expression" << std::endl;
	}

	// Matrix header
	mtx << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	mtx << "%metadata_json: {\"institute\": \"The Crick Institute\"}" << std::endl;
	mtx << GetNumberOfGenes() << " " << GetNumberOfBarcodes() << " " << GetNumberOfEntries() << std::endl;

	// DGE's rows are genes and columns are barcodes. So, we reshape the counts
	struct GeneComp { bool operator() (const Gene& g1, const Gene& g2) const { return g1.GetID() < g2.GetID(); } };
	std::map<Gene, std::map<std::string, unsigned long long>, GeneComp> dge;
	std::map<Gene, unsigned long long>::iterator it;
	for (auto& ent : counts)
	{
		std::map<Gene, unsigned long long> cnts = ent.second.GetCounts();
		for (it = cnts.begin(); it != cnts.end(); ++it)
		{
			dge[it->first][ent.first] += it->second;
		}
	}

	// Matrix
	for (auto& gene : genes_v)
	{
		std::map<std::string, unsigned long long> cnts = dge[gene];
		for (auto& [barcode, cnt] : cnts)
		{
			mtx << positions[gene] << " " << pos[barcode] << " " << cnt << std::endl;
		}
	}
}

// ============================================================================
// Operators
// ============================================================================

// ----------------------------------------------------------------------------
// ostream
// ----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const Counter& counter)
{
	out << "(\n";
	unsigned long long i = 0;
	for (auto& [barcode, counts] : counter.counts)
	{
		out << barcode << ": " << counts;
		i++;
		if ( i < counter.counts.size() )
		{
			out << ",";
		}
		out << "\n";
	}
	out << ")";
	return out;
}

