#!/usr/bin/env python3

import sys
import pandas as pd
from matplotlib.figure import Figure

import matplotlib as mpl
mpl.rc('font', size=16)

###########################
def saturation_plot(x, n):#
###########################
	"""Plots the fraction of reads per cell barcode."""

	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	ax.plot(x.to_list())
	ax.set_ylim((0, 1))
	ax.set_xlabel("Top 10% of cell barcodes, by number of mapped reads")
	ax.set_ylabel("Cumulative fraction of reads")
	
	title = "Cumulative fraction of reads per cell barcode"
	suptitle = "({:,} reads in total)".format(n)
	ax.set_title( title + "\n" + suptitle )

	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path, header=None)
	df.columns = ["Barcode", "UMI", "Reads"]
	reads = df\
		.groupby("Barcode")\
		.apply(lambda x: x.Reads.sum())\
		.to_frame()\
		.reset_index()\
		.rename(columns={0: "Reads"})\
		.sort_values("Reads", ascending=False)
	
	total = reads.Reads.sum()
	
	saturation = reads\
		.iloc[ : reads.shape[0] // 10 ]\
		.Reads\
		.cumsum()\
		.apply(lambda x: x/total)
	
	plt = saturation_plot(saturation, total)
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

