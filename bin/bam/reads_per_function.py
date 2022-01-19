#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#reads_per_function
#samtools view \
#	--expr '[bs] == "MATCHED" && ([cs] == "UNIQUE" || [cs] == "INCLUDED")' \
#	$bam \
#	| tr "\\t" "\\n" \
#	| grep "gf:Z:" \
#	| sed 's/gf:Z://g' \
#	| sort \
#	| uniq -c \
#	| sort -rn \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	> "${name}.reads_per_function.csv"

############################################
def reads_per_function(bam_path, csv_path):#
############################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	functions = []
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		CS = record.get_tag("cs")
		if BS == "MATCHED" and CS in ["UNIQUE", "INCLUDED"]:
			GF = record.get_tag("gf")
			functions.append(GF)
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	df = pd\
		.Series(functions)\
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

