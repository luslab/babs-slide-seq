#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#count_resolved
#samtools view \
#	--expr '[as] == "MAPPED" && [bs] == "MATCHED"' \
#	$bam \
#	| sed 's/^\\([^\\t]\\+\\).*mm:Z:\\([A-Z]\\+\\).*/\\1 \\2/g' \
#	| sort \
#	| uniq \
#	| awk '{ print \$2 }' \
#	| sort \
#	| uniq -c \
#	| sort -rn \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	> "${name}.count_resolved.csv"

########################################
def count_resolved(bam_path, csv_path):#
########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	reads = set()
	counter = 0

	print("BAM iteration", file=sys.stderr)
	for record in bam.fetch(until_eof=True):
		BS = record.get_tag("bs")
		AS = record.get_tag("as")
		if BS == "MATCHED" and AS == "MAPPED":
			MM = record.get_tag("mm")
			reads.add( (record.query_name, MM) )
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)

	print("BAM iteration done, now counting", file=sys.stderr)
	counts = defaultdict(lambda: 0)
	for counter, (r, status) in enumerate(reads):
		counts[status] = counts[status] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	
	print("Export results", file=sys.stderr)
	rows = []
	for k, v in counts.items():
		rows.append({"Status": k, "Reads": v})
	df = pd.DataFrame.from_records(rows)
	df = df.sort_values("Reads", ascending=False)
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	count_resolved(path, csv)

