#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

################################################
def reads_barcode_matching(bam_path, csv_path):#
################################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = defaultdict(set)
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		if BS == "MATCHED":
			barcode = record.get_tag("bm")
			umi = record.get_tag("mi")
			gene = record.get_tag("XF")
			reads[(barcode, umi, gene)].add(record.query_name)
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	counts = defaultdict(lambda:0)
	for k, v in reads.items():
		counts[len(v)] += 1
	
	print("Export results", file=sys.stderr)
	df = pd\
		.Series(counts)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Reads", 0: "Counts"})\
		.sort_values("Counts", ascending=False)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_barcode_matching(path, csv)

