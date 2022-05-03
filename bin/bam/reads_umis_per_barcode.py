#!/usr/bin/python

# coding: utf-8

import sys
import pysam
import pandas as pd
from collections import defaultdict

################################################
def reads_umis_per_barcode(bam_path, csv_path):#
################################################

	bam = pysam.AlignmentFile(bam_path, "rb")
	
	counts = defaultdict(lambda: 0)
	umis = {}
	counter = 0

	
	# TODO: check if at least one barcode passes the UMI threshold
	for record in bam.fetch(until_eof=True):
		up = record.get_tag("us")
		align = record.get_tag("as")
		thresh = record.get_tag("bt")
		if up == "MATCHED" and align == "MAPPED" and thresh == "PASS":
			barcode = record.get_tag("bc")
			counts[barcode] = counts[barcode] + 1
			n_umis = int( record.get_tag("bn") )
			umis[barcode]  = n_umis
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter)
	
	rows = []

	for k in counts.keys():
		rows.append( {"Barcode": k, "Reads": counts[k], "UMIs": umis[k]} )

	df = pd\
		.DataFrame\
		.from_records(rows)\
		.sort_values("UMIs", ascending=False)
	
	df.to_csv(csv_path, header=False, index=False)
	############################################################################

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	csv = sys.argv[2]
	reads_umis_per_barcode(path, csv)

