#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

#reads_per_barcode
#samtools view \
#	--expr '[us] == "MATCHED" && [as] == "MAPPED" && [bc]' $bam \
#	| awk '{
#		v=\$0; gsub(".*bc:Z:", "", v);
#		gsub("\\t.*", "", v);
#		printf "%s\\t%s\\n", v, \$1 }
#	' \
#	| sort \
#	| uniq \
#	| awk '{ print \$1 }' \
#	| sort \
#	| uniq -c \
#	| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
#	| sort -t , -k 2,2rn \
#	> "${csv}"

###########################################
def reads_per_barcode(bam_path, csv_path):#
###########################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	
	counts = defaultdict(lambda: 0)
	counter = 0
	
	for record in bam.fetch(until_eof=True):
		up_status = record.get_tag("us")
		align_status  = record.get_tag("as")
		barcode = record.get_tag("bc")
		if up_status == "MATCHED" and align_status == "MAPPED":
			counts[barcode] = counts[barcode] + 1
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter)
	
	df = pd\
		.Series(counts)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Barcode", 0: "Reads"})\
		.sort_values("Reads", ascending=False)
	
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_per_barcode(path, csv)

