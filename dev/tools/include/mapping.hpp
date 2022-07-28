/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __MAPPING_HPP
#define __MAPPING_HPP

#include <string>
#include <vector>

class Mapping
{
	private:
		std::string gene;
		std::vector<int> scores;
	
	public:
		Mapping();
		Mapping(char*, int);
		Mapping(std::string, int);
		Mapping(const Mapping&);
		std::string GetGene() const;
		std::vector<int> GetScores() const;
		std::string GetScoreString() const;
		int GetScore() const;
		friend std::ostream& operator<<(std::ostream&, const Mapping&);
		friend bool operator<(const Mapping&, const Mapping&);
		friend Mapping operator+(const Mapping&, const Mapping&);
		std::vector<int>::const_iterator begin() const;
		std::vector<int>::const_iterator end() const;
};

#endif

