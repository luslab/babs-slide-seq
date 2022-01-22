#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

##########################################
def reads_per_status(bam_path, csv_path):#
##########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = defaultdict(set)
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		status = record.get_tag("bt")
		reads[status].add(record.query_name)
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	counts = {}
	total = 0
	for k, v in reads.items():
		counts[k] = len(v)
		total = total + counts[k]
	counts["TOTAL"] = total
	
	print("BAM iteration done, now counting", file=sys.stderr)
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
	reads_per_status(path, csv)

