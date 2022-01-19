#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

#########################
def extraction_plot(df):#
#########################
	"""Plots the number of total reads, too short reads and usable reads. The
	argument data frame should have at least 2 columns named "Length" and
	"Reads". The "Length column has to contain an entry named "Total". The
	function returns the modified data frame and a Figure object."""

	# Percentage of total reads
	total = df.Reads.loc[ df.Length == "Total" ].iloc[0]
	df["Percent"] = df.Reads / total * 100

	# annotation
	num = df.Reads.apply(lambda x: "{:,}".format(x)).astype(str)
	percent = df.Percent.round(2).astype(str) + " %"
	df["Annot"] = num + "\n" + percent

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Length,
		"height": df.Reads,
		"color": "skyblue",
		"edgecolor": "black"
	}
	ax.bar(**args)

	# add number and percentage in the middle of the bar
	for i, annot in enumerate(df.Annot):
		args = {
			"s": annot,
			"xy": (i, df.iloc[i].Reads+100),
			"ha": "center",
			"va": "bottom",
			"size": 16,
			"color": "black"
		}
		ax.annotate(**args)

	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/10)
	ax.ticklabel_format(axis="y", style="sci")
	ax.set_title("Barcodes extraction")
	
	fig.tight_layout()

	return df, fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path, header=None, names=["Sample", "Length", "Reads"])
	df = df.filter(["Length", "Reads"])
	df["Length"] = df\
		.Length\
		.replace({
			"Total reads": "Total",
			"Too short reads": "Less than\n41 bp",
			"Long enough reads": "At least\n41 bp"
		})
	df = df\
		.set_index("Length")\
		.loc[ ["Total", "Less than\n41 bp", "At least\n41 bp"] ]\
		.reset_index()
	
	ddf, plt = extraction_plot(df)
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

