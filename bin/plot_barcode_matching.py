#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

###########################
def count_plot(df, title):#
###########################
	"""Plots counts. The argument data frame should have at least 2 columns named
	"Name" and "Count". The function returns a Figure object."""

	# annotation
	df["Annot"] = df.Count.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot.astype(str) + " barcodes"

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Name,
		"height": df.Count,
		"color": "skyblue",
		"edgecolor": "black"
	}
	ax.bar(**args)

	# add number and percentage in the middle of the bar
	for i, annot in enumerate(df.Annot):
		args = {
			"s": annot,
			"xy": (i, df.iloc[i].Count + df.Count.max()/50),
			"ha": "center",
			"va": "bottom",
			"size": 16,
			"color": "black"
		}
		ax.annotate(**args)

	ax.set_ylim(0, df.Count.max() + df.Count.max()/10)
	ax.set_title(title)
	
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path)
	metrics = ["Sequencing barcodes", "Puck barcodes", "Matches"]
	df = df.loc[ df.Metric.isin(metrics) ]
	df.columns = ["Name", "Count"]
	df["Name"] = ["Sequencing", "Puck", "Matched"]
	
	plt = count_plot(df, "Barcode matching")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

