#!/usr/bin/env python3

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

##########################################
def count_duplicates(bam_path, csv_path):#
##########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = {}
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		US = record.get_tag("us")
		AS = record.get_tag("as")
		if US == "MATCHED" and AS == "MAPPED":
			DS = record.get_tag("ds")
			reads[record.query_name] = DS
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	df = pd\
		.Series(reads.values())\
		.value_counts()\
		.reset_index()\
		.rename(columns={"index": "Duplicated", 0: "Reads"})\
		.sort_values("Reads", ascending=False)

	print("Export results", file=sys.stderr)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	count_duplicates(path, csv)

