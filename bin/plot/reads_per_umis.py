#!/usr/bin/python

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

csv_path = "results/qc/anna.reads_per_umi.csv"
out_path = "tmp/hist.pdf"

df = pd.read_csv(csv_path, header=None, names=["Reads", "Elements"])
index = pd.Series(range(1, df.Reads.max()+1))
df = df\
	.set_index("Reads")\
	.reindex(index, fill_value=0)\
	.reset_index()\
	.rename(columns={"index": "Reads"})\
	.assign(Freq=lambda x: x.Elements / x.Elements.sum() * 100)\
	.assign(CumSum=lambda x: x.Freq.cumsum() )

fig = Figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.plot(df.CumSum)
ax.set_ylim((0, 100))
#ax.set_xscale("log")
#ax.set_xlabel(f"Number of {title} (log10)")
#ax.set_title(f"Histogram of {title} per matched barcode")
fig.tight_layout()

fig.savefig(out_path)

###########################
#if __name__ == "__main__":
#
#	mtx_path = sys.argv[1]
#	base_path = sys.argv[2]
#
#	adata = read_10x_mtx(mtx_path)
#	
#	umis = adata\
#		.to_df()\
#		.apply(sum, axis=1)\
#		.to_frame()\
#		.reset_index()\
#		.rename(columns={"index": "Barcode", 0: "UMIs"})
#	
#	THRESH = 10
#	
#	umis = umis.loc[ umis.UMIs >= THRESH ]
#	
#	plt = histo(umis.UMIs, "UMIs")
#	plt.savefig(f"{base_path}.png")
#	plt.savefig(f"{base_path}.pdf")

