#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
from scanpy import read_10x_mtx

import matplotlib as mpl
mpl.rc('font', size=16)

###########################
def count_plot(df, title):#
###########################
	"""Plots counts. The argument data frame should have at least 2 columns named
	"Name" and "Count". The function returns a Figure object."""

	# annotation
	df["Annot"] = df.Count.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot.astype(str) + " UMIs"

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

	mtx_path = sys.argv[1]
	base_path = sys.argv[2]
	
	adata = read_10x_mtx(mtx_path)
	
	umis = adata\
		.to_df()\
		.apply(sum, axis=1)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Barcode", 0: "UMIs"})
	
	top10 = umis.loc[ umis.UMIs >= umis.UMIs.quantile(.9) ]
	
	mean_umis = pd\
		.DataFrame({
			"Name": ["All", "Top 10 %"],
			"Count": [umis.UMIs.mean(), top10.UMIs.mean()]
		})
	
	plt = count_plot(mean_umis, "Number of UMIs per matched barcodes")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

