#!/usr/bin/env python3

# coding: utf-8

import sys, re
import pysam
import pandas as pd
from collections import defaultdict

####################################
def count_tags(bam_path, csv_path):#
####################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	counts = defaultdict(lambda: 0)
	counter = 0
	
	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		XF = record.get_tag("XF")
		if XF.startswith("__"):
			tag = re.sub("\[.*", "", XF)
			counts[tag] = counts[tag] + 1
		else:
			counts["__success"] = counts["__success"] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
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
	count_tags(path, csv)


