#!/usr/bin/python

# coding: utf-8

import sys, re
from os.path import join
import gzip
from Bio import SeqIO

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

###########################
#if __name__ == "__main__":
#	path = sys.argv[1]
#	csv = sys.argv[2]
#	reads_per_barcode(path, csv)

directory = "work/05/8a00703d1a850a26bc7fd5745f034a"
reads1 = [
	"NOO4516A2_S121_L001_R1_001.fastq.gz",
	"NOO4516A2_S121_L002_R1_001.fastq.gz",
	"NOO4516A2_S216_L001_R1_001.fastq.gz",
	"NOO4516A2_S216_L002_R1_001.fastq.gz"
]
reads1 = [
	"NOO4516A2_S121_L001_R2_001.fastq.gz",
	"NOO4516A2_S121_L002_R2_001.fastq.gz",
	"NOO4516A2_S216_L001_R2_001.fastq.gz",
	"NOO4516A2_S216_L002_R2_001.fastq.gz"
]
files = [ join(directory, f) for f in reads1 ]

	
path = "tmp/test.fastq.gz"
fastq = gzip.open(path, "wt")

counter = 0
for f in files:
	handle = gzip.open(f, "rt")
	fq = SeqIO.parse(handle, "fastq")
	for read in fq:
		fastq.write(read.format("fastq"))
		counter = counter + 1
		if counter % 1000000 == 0:
			print(counter, file=sys.stderr)
	handle.close()
fastq.close()

