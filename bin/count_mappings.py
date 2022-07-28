#!/usr/bin/env python3

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

########################################
def count_mappings(bam_path, csv_path):#
########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	genes = defaultdict(set)
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		genes[record.query_name].add( record.get_tag("qi") )
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	counts = {}
	for k, v in genes.items():
		counts[k] = len(v)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	df = pd\
		.Series(counts.values())\
		.value_counts()\
		.reset_index()\
		.rename(columns={"index": "Mapping", 0: "Reads"})\
		.sort_values("Reads", ascending=False)

	print("Export results", file=sys.stderr)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	count_mappings(path, csv)

