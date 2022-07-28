#!/usr/bin/env python3

import sys
from random import randint
from io import StringIO
import pandas as pd
import gzip
from Bio import SeqIO

###################################################################################
def duplicates(path, fastq_1, fastq_2, csv, molecule_delim="|", record_delim="/"):#
###################################################################################

	n = sum(1 for line in open(path, "r"))
	nums = [ randint(0, n) for i in range(20) ]
	rows = [ line.rstrip() for i, line in enumerate(open(path, "r")) if i in nums ]
	
	df = pd\
		.read_table(StringIO("\n".join(rows)), sep=",", names=["Barcode", "UMI", "Reads"])\
		.assign(Reads=lambda x: x.Reads.map(lambda x: x.split(molecule_delim)))\
		.explode("Reads")\
		.assign(Reads=lambda x: x.Reads.map(lambda x: x.split(record_delim)))\
		.assign(Read=lambda x: x.Reads.map(lambda x: x[0]))\
		.assign(Gene=lambda x: x.Reads.map(lambda x: x[1]))\
		.assign(Score=lambda x: x.Reads.map(lambda x: x[2]))\
		.drop("Reads", axis=1)\
		.assign(Read1="", Read2="")\
		.set_index("Read")
	
	fq = {1: fastq_1, 2: fastq_2}
	for k, v in fq.items():
		f = gzip.open(v, "rt")
		i = 0
		for record in SeqIO.parse(f, "fastq"):
			if record.id in df.index:
				df.loc[record.id]["Read"+str(k)] = str(record.seq)
			i += 1
			if i % 1000000 == 0:
				print(i)
		f.close()
	
	df = df.reset_index()
	df = df.sort_values(["Barcode", "UMI"])
	df.to_csv(csv, index=False)
	###########################################################################

#args = {
#	"molecule_delim": "|",
#	"record_delim": "/",
#	"path": "results/files/anna.select.resolved.csv",
#	"fastq_1": "results/files/anna_R1.fastq.gz",
#	"fastq_2": "results/files/anna_R2.fastq.gz",
#	"csv": "tmp/anna_dup.csv"
#	}

##########################
if __name__ == "__main__":
	path = sys.argv[1]
	fastq_1 = sys.argv[2]
	fastq_2 = sys.argv[3]
	csv = sys.argv[4]
	duplicates(path, fastq_1, fastq_2, csv)

