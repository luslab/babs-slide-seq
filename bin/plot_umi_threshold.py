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
	"""Plots counts. The argument data frame should have at least 2 columns
	named "Name" and "Count". The function returns a Figure object."""

	# annotation
	df["Percent"] = np.round(df.Reads / df.Reads.sum() * 100, 1)
	df["Annot"] = df.Reads.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot + " reads\n" + df.Percent.astype(str) + " %"

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Status,
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
			"size": 16,
			"color": "black"
		}
		ax.annotate(**args)

	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/8)
	suptitle = "({:,} reads in total)".format(df.Reads.sum())
	ax.set_title( title + "\n" + suptitle )
	
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	threshold = int( sys.argv[3] )

	df = pd.read_csv(csv_path, header=None)
	df.columns = ["Status", "Reads"]
	status = {
		"PASS": f"$\geq${threshold} UMIs",
		"TOO_LOW": f"$<${threshold} UMIs"
	}
	df = df.loc[ df.Status.isin( status.keys() ) ]
	df = df\
		.set_index("Status")\
		.reindex(status.keys(), fill_value=0)\
		.reset_index()
	df["Status"] = df.Status.map(status)
	df["Status"] = pd.Categorical(df.Status, categories=status.values())
	
	plt = count_plot(df, "Reads with barcode that pass the UMIs threshold")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

