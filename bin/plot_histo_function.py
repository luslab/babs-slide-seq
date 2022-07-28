#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

##############################
def function_plot(df, title):#
##############################
	"""Plots the number reads per function. The argument data frame should have
	at least 2 columns named "Function" and "Reads". The function returns the
	modified data frame and a Figure object."""

	coding = df.loc[df.Function == "CODING"].iloc[0]["Reads"]
	utr = df.loc[df.Function == "UTR"].iloc[0]["Reads"]
	total = df.Reads.sum()

	# Percentage of total reads
	df["Percent"] = df.Reads / total * 100

	# annotation
	df["Percent"] = np.round(df.Reads / df.Reads.sum() * 100, 1)
	df["Annot"] = df.Reads.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot + " reads\n" + df.Percent.astype(str) + " %"

	df["Function"] = df.Function.str.title()
	df["Function"] = df.Function.str.replace("Utr", "UTR")

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Function,
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
	ax.set_xticklabels(df.Function, size=14)

	sub_title = "({:,} Coding and UTR/{:,} reads)".format(coding+utr, total)
	ax.set_title( title + "\n" + sub_title , size=18 )
	
	fig.tight_layout()

	return df, fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path, header=None, names=["Function", "Reads"])
	df["Function"] = df.Function.replace({np.nan: "Ambiguous"})
	
	ddf, plt = function_plot(df, "Gene functions for bead-matched reads")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

