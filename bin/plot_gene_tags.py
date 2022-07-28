#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

#########################
def tag_plot(df, title):#
#########################
	"""Plots the number reads per tag. The argument data frame should have at
	least 2 columns named "Tag" and "Reads". The function returns the modified
	data frame and a Figure object."""

	success = df.loc[df.Tag == "Success"].iloc[0]["Reads"]
	total = df.Reads.sum()

	# Percentage of total reads
	df["Percent"] = df.Reads / total * 100

	# annotation
	df["Percent"] = np.round(df.Reads / df.Reads.sum() * 100, 1)
	df["Annot"] = df.Reads.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot + " reads\n" + df.Percent.astype(str) + " %"

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Tag,
		"height": df.Reads,
		"color": "skyblue",
		"edgecolor": "black"
	}
	ax.bar(**args)

	# add number and percentage in the middle of the bar
	for i, annot in enumerate(df.Annot):
		args = {
			"s": annot,
			"xy": (i, df.iloc[i].Reads + df.Reads.max()/50),
			"ha": "center",
			"va": "bottom",
			"size": 12,
			"color": "black"
		}
		ax.annotate(**args)

	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/10)
	ax.ticklabel_format(axis="y", style="sci", scilimits=(6,6))
	ax.set_xticks(range(df.Reads.size))
	ax.set_xticklabels(df.Tag, size=12)

	sub_title = "({:,} good/{:,} reads)".format(success, total)
	ax.set_title( title + "\n" + sub_title , size=18 )
	
	fig.tight_layout()

	return df, fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	names = {
		"__success": "Success",
		"__no_feature": "No feature",
		"__ambiguous": "Ambiguous",
		"__too_low_aQual": "Low quality",
		"__not_aligned": "Not Aligned",
		"__alignment_not_unique": "Not unique"
		}
	
	df = pd.read_csv(csv_path, header=None, names=["Tag", "Reads"])
	df["Tag"] = df.Tag.map(names)
	
	ddf, plt = tag_plot(df, "Gene annotation of the reads")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

