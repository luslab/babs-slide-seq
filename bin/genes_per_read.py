#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

########################################
def genes_per_read(bam_path, csv_path):#
########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = {}
	counter = 0
	
	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		reads[record.query_name] = record.get_tag("rm")
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("Counting")
	df = pd\
		.Series(reads)\
		.value_counts()\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Status", 0: "Reads"})
	
	print("Exporting")
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	genes_per_read(path, csv)

