#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

###########################
def genes_plot(df, title):#
###########################
	"""Plots number of matched reads. The argument data frame should have at
	least 2 columns named "Status" and "Reads"."""

	one_gene = df.loc[df.Status=="One gene"].iloc[0]["Reads"]
	several_genes = df.loc[df.Status=="Several genes"].iloc[0]["Reads"]
	with_gene = one_gene + several_genes

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
	suptitle = "({:,} with genes/{:,} reads)".format(with_gene, df.Reads.sum())
	ax.set_title( title + "\n" + suptitle )
	
	fig.tight_layout()

	return fig
	############################################################################
	
##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]

	df = pd.read_csv(csv_path, header=None, names=["Status", "Reads"])
	status = {
		"ONE_GENE": "One gene",
		"SEVERAL_GENES": "Several genes",
		"UNASSIGNED": "Unassigned"
	}
	df = df.loc[ df.Status.isin( status.keys() ) ]
	df = df\
		.set_index("Status")\
		.reindex(status.keys(), fill_value=0)\
		.reset_index()
	df["Status"] = df.Status.map(status)
	df["Status"] = pd.Categorical(df.Status, categories=status.values())
		
	plt = genes_plot(df, "Gene mapping for reads")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

