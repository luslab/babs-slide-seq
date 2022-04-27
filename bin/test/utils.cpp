/*
 * nourdine.bah@crick.ac.uk
 */

#include <iostream>
#include <map>
#include <string>

#include "utils.hpp"

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	////////////////////////////////////////////////////////////////////////////
	std::string seq1 = "AAAGC";
	std::string seq2 = "AAATC";
	std::cout << seq1 << ", " << seq2 << ": " << HammingDist(seq1, seq2) << std::endl;

	////////////////////////////////////////////////////////////////////////////
	std::cout << DummySeq(8) << std::endl;

	////////////////////////////////////////////////////////////////////////////
	std::string sample = "Sample";

	std::map<std::string, unsigned long long> metrics;
	metrics["Metric1"] = 6;
	metrics["Metric2"] = 12;
	std::string metrics_csv = "/tmp/metrics.csv";
	WriteCounter<std::string, unsigned long long>(sample, metrics_csv, metrics);

	std::map<int, float> histo;
	histo[5] = 4.5;
	histo[100] = 89.45;
	std::string histo_csv = "/tmp/histo.csv";
	WriteCounter<int, float>(sample, histo_csv, histo);

	return 0;
}

