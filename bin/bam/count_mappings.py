#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#count_mappings
#samtools view \
#	--expr '[bs] == "MATCHED" && [as] == "MAPPED" && [ds] == "PRIMARY"' \
#	$bam \
#	| awk '{ print \$1 }' \
#	| sort \
#	| uniq -c \
#	| awk '{ print \$1 }' \
#	| sort \
#	| uniq -c \
#	| sort -rn \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	> "${name}.count_mappings.csv"

########################################
def count_mappings(bam_path, csv_path):#
########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	counts = defaultdict(lambda: 0)
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		AS = record.get_tag("as")
		DS = record.get_tag("ds")
		if BS == "MATCHED" and AS == "MAPPED" and DS == "PRIMARY":
			counts[record.query_name] = counts[record.query_name] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
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

