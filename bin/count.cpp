#include <map>
#include <set>
#include <seqan/basic.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>
#include <seqan/store.h>

using namespace seqan;

// for gff
typedef std::tuple<CharString, CharString, int> Gene;
typedef std::set<Gene> GeneSet;
typedef std::vector<Gene> Genes;
typedef std::map<CharString, Gene> GeneMap;

// for bam
typedef std::set<CharString> Umis;
typedef std::map<CharString, Umis> Barcodes;
typedef std::map<CharString, Barcodes> GeneCount;

///////////////////////////////////////////////////////////////////////////////

//////////////////////////////////
bool compareGene(Gene g1, Gene g2)
{
	CharString id1 = std::get<1>(g1);
	CharString id2 = std::get<1>(g2);
	std::string s1(toCString(id1));
	std::string s2(toCString(id2));
	return (s1<s2);
}

////////////////////////////////////////////////////////////////////
int get_position(StringSet<CharString> tagNames, CharString tagName)
{
	int position = -1;

	for (int i=0; i<length(tagNames); i++)
	{
		if ( tagName == tagNames[i] )
		{
			position = i;
		}
	}

	return position;
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	std::string out_dir(argv[3]);

	////////////////////////////////////////////////////////////////////////////
	GffFileIn gff(argv[1]);
	GffRecord record;

	GeneSet gene_set;

	////////////////////////////////////////////////////////////////////////////
	std::set<CharString> TYPES = {
		"gene", "exon",
		"start_codon", "stop_codon",
		"five_prime_utr", "three_prime_utr"
	};

	std::set<CharString>::iterator t_it;

	////////////////////////////////////////////////////////////////////////////
	
	int g = 1;

	std::map<CharString, CharString> names;

	while ( ! atEnd(gff) )
	{
		readRecord(record, gff);

		/////////////////////////////////////////////////////////////////////////
		t_it = TYPES.find(record.type);
		if ( t_it == TYPES.end() )
		{
			continue;
		}
		/////////////////////////////////////////////////////////////////////////

		int namePos = get_position(record.tagNames, "gene_name");
		int idPos = get_position(record.tagNames, "gene_id");

		if ( -1 != namePos && -1 != idPos )
		{
			Gene gene = std::make_tuple(
				record.tagValues[idPos],
				record.tagValues[namePos],
				0
			);
			gene_set.insert(gene);

			names[ record.tagValues[idPos] ] = record.tagValues[namePos];
		}
		// counter
		if ( g % (unsigned long)100000 == 0 ) {
			std::cerr << g << " gff entries" << std::endl;
		}
		g++;
	}

	Genes genes(gene_set.begin(), gene_set.end());
	std::sort(genes.begin(), genes.end(), compareGene);
	GeneMap gene_map;
	for (int i=0; i<genes.size(); i++)
	{
		Gene gene = std::make_tuple(
			std::get<0>(genes[i]),
			std::get<1>(genes[i]),
			i+1
		);

		gene_map[ std::get<0>(genes[i]) ] = gene;
	}

	////////////////////////////////////////////////////////////////////////////

	BamFileIn bam(argv[2]);
	BamHeader header;
	readHeader(header, bam);
	BamAlignmentRecord rec;

	CharString barcode_status;
	CharString counting_status;

	CharString barcode;
	CharString umi;
	CharString symbol;
	CharString id;

	unsigned int tag_idx = 0;

	unsigned long counter = 0;

	std::set<CharString> barcode_set;
	GeneCount count;
	unsigned long n_entry = 0;

	while ( ! atEnd(bam) )
	{
		readRecord(rec, bam);
		BamTagsDict tagsDict(rec.tags);

		if ( findTagKey(tag_idx, tagsDict, "bc") ) {
			extractTagValue(barcode, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find barcode" << std::endl;
			return 1;
		}

		if ( findTagKey(tag_idx, tagsDict, "mi") ) {
			extractTagValue(umi, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find UMI" << std::endl;
			return 1;
		}

		if ( findTagKey(tag_idx, tagsDict, "XF") ) {
			extractTagValue(id, tagsDict, tag_idx);
		}
		else
		{
			std::cerr << "Problem: cannot find gene ID" << std::endl;
			return 1;
		}

		if ( count.find(id) == count.end() )
		{
			Umis umis;
			Barcodes barcodes;
			umis.insert(umi);
			barcodes[barcode] = umis;
			count[id] = barcodes;

			n_entry++;
		}
		else
		{
			if ( count[id].find(barcode) == count[id].end() )
			{
				Umis umis;
				umis.insert(umi);
				count[id][barcode] = umis;

				n_entry++;
			}
			else
			{
				count[id][barcode].insert(umi);
			}
		}

		barcode_set.insert(barcode);

		/*
		std::cout
			<< status
			<< ", "
			<< barcode
			<< ", "
			<< umi
			<< ", "
			<< id
			<< ", "
			<< symbol
			<< std::endl;
		// */

		/////////////////////////////////////////////////////////////////////////
		// counter
		counter++;
		if ( counter % (unsigned long)1000000 == 0 ) {
			std::cerr << counter << std::endl;
		}
		//if ( counter > 1000000 ) { break; }
	}

	close(bam);

	////////////////////////////////////////////////////////////////////////////
	// BEAD ORDERING
	
	std::vector<CharString> barcodes(barcode_set.begin(), barcode_set.end());
	std::sort(barcodes.begin(), barcodes.end());
	std::map<CharString, long unsigned int> barcode_order;

	std::ofstream bcd(out_dir + "/barcodes.tsv");

	for (int i=0; i<barcodes.size(); i++)
	{
		barcode_order[ barcodes[i] ] = i + 1;
		bcd << barcodes[i] << std::endl;
	}

	bcd.close();

	////////////////////////////////////////////////////////////////////////////
	// MATRIX

	std::ofstream feat(out_dir + "/features.tsv");

	std::ofstream mtx(out_dir + "/matrix.mtx");

	mtx << "%%MatrixMarket matrix coordinate integer general" << std::endl;

	mtx
		<< "%metadata_json: {\"institute\": \"The Crick Institute\"}"
		<< std::endl;

	mtx
		<< gene_map.size()
		<< " "
		<< barcodes.size()
		<< " "
		<< n_entry
		<< std::endl;

	for (int i=0; i<genes.size(); i++)
	{
		CharString gene_id = std::get<0>(genes[i]);

		// features
		feat
			<< gene_id
			<< "\t"
			<< std::get<1>(genes[i])
			<< "\t"
			<< "Gene Expression"
			<< std::endl;

		// counting
		if ( count.find(gene_id) != count.end() )
		{
			Barcodes bc = count[gene_id];

			for (int j=0; j<barcodes.size(); j++)
			{
				if ( bc.find( barcodes[j] ) != bc.end() )
				{
					mtx
						<< i+1
						<< " "
						<< j+1
						<< " "
						<< bc[ barcodes[j] ].size()
						<< std::endl;
				}
			}
		}
	}
	
	feat.close();
	mtx.close();

	////////////////////////////////////////////////////////////////////////////
	return 0;
}
