#!/usr/bin/python

# coding: utf-8

import sys, re
import pandas as pd
import gzip
from Bio import SeqIO
from collections import defaultdict

#reads_per_barcode
#zcat $fastq \
#	| sed -n '1~4p' \
#	| grep LONG_ENOUGH \
#	| grep ":MATCHED:" \
#	| sed 's/^[^:]\\+:[^:]\\+:[^:]\\+:[^:]\\+:\\([^:]\\+\\).*/\\1/g' \
#	| sort \
#	| uniq -c \
#	| sort -rh \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	> "${name}.reads_per_barcode.csv"

#############################################
def reads_per_barcode(fastq_path, csv_path):#
#############################################

	counts = defaultdict(lambda: 0)
	counter = 0
	regex = re.compile(r"([A-Z_]+):[A-Z]+:\d+:([A-Z]+):([A-Z]+):.*")
		
	with gzip.open(fastq_path, "rt") as fastq:
		for record in SeqIO.parse(fastq, "fastq"):
			m = regex.match(record.id)
			if m:
				length_status = m.group(1)
				match_status = m.group(2)
				barcode = m.group(3)
				if length_status == "LONG_ENOUGH" and match_status == "MATCHED":
					counts[barcode] = counts[barcode] + 1
			counter = counter + 1
			if counter % 1000000 == 0:
				print(counter, file=sys.stderr)
	
	df = pd.DataFrame({"Barcode": counts.keys(), "Reads": counts.values()})
	df = df.sort_values("Reads", ascending=False)
	df.to_csv(csv_path, index=False, header=False)
	###########################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_per_barcode(path, csv)

