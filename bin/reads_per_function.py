#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

############################################
def reads_per_function(bam_path, csv_path):#
############################################

	order = {
		"CODING": 0,
		"UTR": 1,
		"INTRONIC": 2,
		"INTERGENIC": 3
	}
	
	bam = pysam.AlignmentFile(bam_path, "rb")
	functions = defaultdict(set)
	counter = 0
	
	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		GF = record.get_tag("qf")
		functions[record.query_name].add(GF)
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now classifying", file=sys.stderr)
	classif = {}
	counter = 0
	for k, v in functions.items():
		classif[k] = sorted(list(v), key=lambda x: order[x])[0]
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("Counting", file=sys.stderr)
	df = pd\
		.Series(classif)\
		.value_counts()\
		.reset_index()\
		.rename(columns={"index": "Function", 0: "Reads"})\
		.sort_values("Reads", ascending=False)
	
	print("Export results", file=sys.stderr)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_per_function(path, csv)

