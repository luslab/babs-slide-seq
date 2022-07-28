#!/usr/bin/env python3

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

###########################################
def reads_up_matching(bam_path, csv_path):#
###########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = set()
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		LS = record.get_tag("ls")
		if LS == "LONG_ENOUGH":
			US = record.get_tag("us")
			AS = record.get_tag("as")
			reads.add( ( record.query_name , US , AS ) )
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	counts = defaultdict(lambda: 0)
	for counter, r in enumerate(reads):
		counts[(r[1], r[2])] = counts[(r[1], r[2])] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("Export results", file=sys.stderr)
	rows = []
	for k, v in counts.items():
		rows.append({"Matched": k[0], "Mapped": k[1], "Reads": v})
	df = pd.DataFrame.from_records(rows)
	df = df.sort_values("Reads", ascending=False)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_up_matching(path, csv)

