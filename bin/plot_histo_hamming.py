#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

######################
def hamming_plot(df):#
######################
	"""Plots histgram of hamming distance for shuffled and not shuffled barcodes.
	The function returns a Figure object."""

	fig = Figure(figsize=(8,8))
	ax = fig.add_subplot(111)

	color = {"Ordered": "blue", "Shuffled": "red"}
	
	for barcodes, gdf in df.groupby("Barcodes"):
		ax.bar(
			gdf.Distance, gdf.Count,
			width=1,
			color=color[barcodes],
			edgecolor="black",
			label=barcodes,
			alpha=.3
			)
	
	ax.set_xlabel("Hamming distance")
	ax.set_ylabel("Number of barcodes")
	ax.set_title("Histogram of hamming distances")
	
	ax.legend()
	
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	ordered_path = sys.argv[1]
	shuffled_path = sys.argv[2]
	base_path = sys.argv[3]
	
	columns = {0: "PuckBarcode", 1: "Distance", 2: "Matches", 3: "SeqBarcodes"}
	
	ordered = pd\
		.read_csv(ordered_path, header=None)\
		.rename(columns=columns)\
		.Distance\
		.value_counts()\
		.to_frame()\
		.rename(columns={"Distance": "Count"})
	
	shuffled = pd\
		.read_csv(shuffled_path, header=None)\
		.rename(columns=columns)\
		.Distance\
		.value_counts()\
		.to_frame()\
		.rename(columns={"Distance": "Count"})
	
	max_dist = max(ordered.index.max(), shuffled.index.max()) + 3
	
	ordered = ordered\
		.reindex(range(max_dist), fill_value=0)\
		.assign(Barcodes="Ordered")
	
	shuffled = shuffled\
		.reindex(range(max_dist), fill_value=0)\
		.assign(Barcodes="Shuffled")
	
	hamming = pd\
		.concat([ordered, shuffled])\
		.reset_index()\
		.rename(columns={"index": "Distance"})
	
	plt = hamming_plot(hamming)
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

