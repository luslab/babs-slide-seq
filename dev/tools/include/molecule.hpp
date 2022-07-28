/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __MOLECULE_HPP
#define __MOLECULE_HPP

#include <string>
#include <set>
#include <map>

#include "record.hpp"
#include "mappings.hpp"

class Molecule
{
	private:
		std::string barcode;
		std::string umi;
		std::string sequence;
		std::set<Record> records;
		Mappings mappings;
		std::map<std::string, long long> frequencies;
		std::string status;
	
	public:
		Molecule();
		Molecule(char*, char*, Record);
		Molecule(const Molecule&);
		std::string GetBarcode() const;
		std::string GetUMI() const;
		std::string GetSeq() const;
		std::string GetStatus() const;
		std::set<Record> GetRecords() const;
		std::string GetRecordString() const;
		std::set<std::string> GetGenes() const;
		std::string GetGeneString() const;
		void Insert(Record);
		void ExtractMappings();
		Mappings GetMappings(bool);
		long GetMaxScore() const;
		bool IsThereAMaxima() const;
		std::map<unsigned long long, std::string> GetRecordTags() const;
		void ComputeFrequencies();
		long long GetMaxFrequency() const;
		bool IsThereAMajority() const;
		std::map<unsigned long long, std::string> GetFrequencyBasedRecordTags();
		std::set<std::string> GetReads() const;
		std::string GetCSVString(char, char) const;
		friend std::ostream& operator<<(std::ostream&, const Molecule&);
		friend bool operator<(const Molecule&, const Molecule&);
		friend Molecule operator+(Molecule&, Molecule&);
		std::set<Record>::iterator begin();
		std::set<Record>::iterator end();
};

#endif

