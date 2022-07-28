#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

######################################
def count_select(bam_path, csv_path):#
######################################

	order = {
		"UNIQUE": 0,
		"INCLUDED": 1,
		"EXCLUDED": 2,
		"UNRESOLVED": 3
	}

	bam = pysam.AlignmentFile(bam_path, "rb")
	status = defaultdict(set)
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		status[record.query_name].add( record.get_tag("cs") )
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	counts = defaultdict(lambda: 0)
	for k, v in status.items():
		tag = sorted(list(v), key=lambda x: order[x])[0]
		counts[tag] = counts[tag] + 1
	df = pd\
		.Series(counts)\
		.reset_index()\
		.rename(columns={"index": "Status", 0: "Reads"})\
		.sort_values("Reads", ascending=False)

	print("Export results", file=sys.stderr)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	count_select(path, csv)

