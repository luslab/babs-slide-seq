#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#reads_per_barcode_umi
#samtools view \
#	--expr '[bs] == "MATCHED" && ([cs] == "UNIQUE" || [cs] == "INCLUDED")' \
#	$bam \
#	| sed 's/^\\([^\\t]\\+\\).*bc:Z:\\([A-Z]\\+\\).*mi:Z:\\([A-Z]\\+\\).*/\\2 \\3 \\1/g' \
#	| sort \
#	| cut -d " " -f 1,2 \
#	| uniq -c \
#	| sort -nr \
#	| awk '{ printf "%s,%s,%s\\n", \$2, \$3, \$1 }' \
#	> "${name}.reads_per_barcode_umi.csv"

###############################################
def reads_per_barcode_umi(bam_path, csv_path):#
###############################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = set()
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		CS = record.get_tag("cs")
		if BS == "MATCHED" and CS in ["UNIQUE", "INCLUDED"]:
			BC = record.get_tag("bc")
			MI = record.get_tag("mi")
			reads.add( ( record.query_name , BC , MI ) )
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("BAM iteration done, now counting", file=sys.stderr)
	counts = defaultdict(lambda: 0)
	for counter, (r, barcode, umi) in enumerate(reads):
		counts[(barcode, umi)] = counts[(barcode, umi)] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("Export results", file=sys.stderr)
	rows = []
	for (barcode, umi), cnt in counts.items():
		rows.append({"Barcode": barcode, "UMI": umi, "Reads": cnt})
	df = pd.DataFrame.from_records(rows)
	df = df.sort_values("Reads", ascending=False)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_per_barcode_umi(path, csv)

