/*
 * nourdine.bah@crick.ac.uk
 */

#include <stdio.h>
#include <stdlib.h>

#include <set>
#include <string>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>

#include "data.hpp"

int main(int argc, char** argv)
{
	std::string idx_path("/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp");
	std::string gff_path("/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf");
	std::string fasta_path("/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome/Mus_musculus.GRCm38.dna_sm.toplevel.fa");
	std::string fasta_index_path("/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome/Mus_musculus.GRCm38.dna_sm.toplevel.fa.fai");
	std::string read_structure("8C18U7C2X8M");

	std::string coords_path("fake_puck/slideseq.csv");

	unsigned long long n_beads = 80000;

	// design file
	std::ofstream design("design.csv");
	design << "name,fastq_1,fastq_2,puck,read_structure,genome,gtf" << std::endl;

	for (int sample=1; sample <3; sample++)
	{
		std::cerr << "Creating sample " << sample << "..." << std::endl;
		Data data = Data((char*)fasta_path.c_str(), (char*)fasta_index_path.c_str(),
			(char*)gff_path.c_str(), (char*)coords_path.c_str(), n_beads);


		// barcode coordinates
		std::cerr << "Writing coordinages for sample " << sample << "..." << std::endl;
		std::vector<std::tuple<std::string, double, double>> barcodes = data.GetBarcodes();
		std::ofstream puck("puck" + std::to_string(sample) + ".csv");
		puck << "Barcode,x,y" << std::endl;
		for (unsigned long long b=0; b<barcodes.size(); b++)
		{
			puck << std::get<0>(barcodes[b]) << "," << std::get<1>(barcodes[b]) << "," << std::get<2>(barcodes[b]) << std::endl;
		}
		puck.close();

		for (int lane=1; lane<6; lane++)
		{
			std::cerr << "Creating lane " << lane << "..." << std::endl;

			std::string path1 = "sample" + std::to_string(sample) + "_L00" + std::to_string(lane) + "_R1.fastq.gz";
			seqan::SeqFileOut fastq1(path1.c_str());
			std::string path2 = "sample" + std::to_string(sample) + "_L00" + std::to_string(lane) + "_R2.fastq.gz";
			seqan::SeqFileOut fastq2(path2.c_str());

			data.GetSequences(fastq1, fastq2, 20000, "sample" + std::to_string(sample) + "_lane" + std::to_string(lane) + "_");

			seqan::close(fastq1);
			seqan::close(fastq2);

			design << "sample" + std::to_string(sample) + ",test/" + path1 + ",test/" + path2 + ",test/puck" + std::to_string(sample) + ".csv," + read_structure + "," + idx_path + "," + gff_path << std::endl;
		}
	}

	design.close();

	return 0;
}

