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
def count_plot(df, title):#
###########################
	"""Plots counts. The argument data frame should have at least 2 columns
	named "Name" and "Count". The function returns a Figure object."""

	# annotation
	df["Percent"] = np.round(df.Count / df.Count.sum() * 100, 1)
	df["Annot"] = df.Count.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot + " reads (" + df.Percent.astype(str) + " %)"

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
	ax.set_title( title + "\n" + "({:,} reads)".format(df.Count.sum()) )
	
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]

	df = pd.read_csv(csv_path, header=None)
	df.columns = ["Mappings", "Reads"]
	df = df.set_index("Mappings")
	df = df.reindex(range(1, df.index.max()+1), fill_value=0).reset_index()
	df["Name"] = np.where(df.Mappings > 1, "Multi-mapped", "Uniquely mapped")
	cats = ["Uniquely mapped", "Multi-mapped"]
	df["Name"] = pd.Categorical(df.Name, categories=cats)
	
	mapping = df\
		.groupby("Name")\
		.apply(lambda x: x.Reads.sum())\
		.to_frame()\
		.reset_index()\
		.rename(columns={0: "Count"})
	
	plt = count_plot(mapping, "Multi-mapping for reads")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

