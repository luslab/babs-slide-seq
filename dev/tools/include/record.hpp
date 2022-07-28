/*
 * nourdine.bah@crick.ac.uk
 */

#ifndef __RECORD_HPP
#define __RECORD_HPP

#include <string>

class Record
{
	private:
		unsigned long long position;
		std::string read;
		long score;
		std::string gene;
	
	public:
		Record();
		Record(unsigned long long, char*, long, char*);
		unsigned long long GetPos() const;
		std::string GetRead() const;
		long GetScore() const;
		std::string GetGene() const;
		std::string GetCSVString(char) const;
		friend std::ostream& operator<<(std::ostream&, const Record&);
		friend bool operator<(const Record&, const Record&);
};

#endif

