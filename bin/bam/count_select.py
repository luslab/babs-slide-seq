#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#count_select
#samtools view \
#	--expr '[bs] == "MATCHED" && ( [mm] == "UNIQUE" || [mm] == "INCLUDED" )' \
#	$bam \
#	| sed 's/^.*cs:Z:\\([A-Z]\\+\\).*/\\1/g' \
#	| sort \
#	| uniq -c \
#	| sort -rn \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	> "${name}.count_select.csv"

######################################
def count_select(bam_path, csv_path):#
######################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	status = []
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		MM = record.get_tag("mm")
		if BS == "MATCHED" and MM in ["UNIQUE", "INCLUDED"]:
			CS = record.get_tag("cs")
			status.append(CS)
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	df = pd\
		.Series(status)\
		.value_counts()\
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

