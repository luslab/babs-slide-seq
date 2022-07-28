#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scanpy import read_10x_mtx

import matplotlib as mpl
mpl.rc('font', size=16)

#####################
def histo(x, title):#
#####################
	"""Plots an histogram."""

	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	ax.hist(
		x,
		bins=np.logspace(0, 4, 41),
		facecolor="skyblue",
		edgecolor="black"
	)
	ax.set_xscale("log")
	ax.set_xlabel(f"Number of {title} (log10)")
	ax.set_title(f"Histogram of {title} per matched barcode")
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	mtx_path = sys.argv[1]
	base_path = sys.argv[2]

	adata = read_10x_mtx(mtx_path)
	
	genes = adata\
		.to_df()\
		.apply(lambda x: x.apply(lambda y: 1 if y>0 else 0).sum(), axis=1)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Barcode", 0: "Genes"})
	
	THRESH = 10
	
	genes = genes.loc[ genes.Genes >= THRESH ]
	
	plt = histo(genes.Genes, "Genes")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

