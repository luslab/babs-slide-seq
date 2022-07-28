#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

#############################################
def base_balance_plot(df, title, sub_title):#
#############################################
	"""Plots the base balance. The argument data frame should have its columns as
	ligations and bases as index. The function returns a Figure object."""

	colors = ["red", "green", "blue", "purple", "grey"]

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)

	for i, base in enumerate(df.index):

		args = {
			"x": df.loc[base].index,
			"y": df.loc[base],
			"s": 22,
			"color": colors[i],
			"label": base
		}

		ax.scatter(**args)
	
	ax.set_xlabel("Ligation")
	ax.set_ylabel("% of reads")
	ax.legend(prop={"size":12}, ncol=2)

	ax.set_xticks(df.loc[base].index)
	ax.set_xticklabels(df.loc[base].index, size=14)

	ax.set_title( title + sub_title , size=18 )
	
	fig.tight_layout()

	return fig
	############################################################################

#####################
def base_balance(x):#
#####################
	balance = {}
	for base in ["A", "C", "G", "T", "N"]:
		balance[base] = np.sum( x == base ) / x.size * 100
	return pd.Series(balance)
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd\
		.read_csv(csv_path, header=None)\
		.rename(columns={0: "Barcode", 1: "UMI", 2: "Reads"})
	
	barcodes = pd.Series(df.Barcode.sort_values().unique())
	barcodes = pd.DataFrame(barcodes.apply(list).to_list())
	barcodes.columns = [ c + 1 for c in barcodes.columns ]
	barcodes = barcodes.apply(base_balance, axis=0)
	
	plt = base_balance_plot(
		barcodes,
		"Barcode for bead-matched reads",
		"\n({:,} reads)".format(df.Reads.sum())
		)
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

