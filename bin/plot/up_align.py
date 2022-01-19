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
def align_plot(df, title):#
###########################
	"""Plots number of mapped and matched reads. The argument data frame should
	have at least 3 columns named "Mapped", "Matched" and "Reads. The "Mapped"
	and "Matched columnms has to be of binary values."""

	matched = df.loc[ df.Matched ].Reads.sum()
	total = df.Reads.sum()

	# percentage and ordering
	df["Percent"] = df.Reads / total * 100
	df = df.sort_values(["Mapped", "Matched"])
	df.index = range(df.shape[0])

	# the plot
	fig = Figure(figsize=(8,8))
	ax = fig.add_subplot(111)
	width = .4
	xticklabels = []
	
	for i, row in df.iterrows():
	
		x = 1 if row.Mapped else 0
		pos = 1 if row.Matched else -1
		color = "tab:orange" if row.Matched else "skyblue"
		label = "Matched" if row.Matched else "Unmatched"
		xlabel = "Mapped" if row.Mapped else "Unmapped"
	
		args = {
			"x": x + pos * width / 2,
			"height": row.Reads,
			"width": width,
			"color": color,
			"edgecolor": "black"
		}
	
		if row.Mapped:
			args = {**args, **{"label": label}}
		
		ax.bar(**args)
	
		s = "{:,}".format(row.Reads) + "\n" + str(round(row.Percent, 2)) + " %"
		args = {
			"s": s,
			"xy": (x + pos * width / 2, row.Reads + 100),
			"ha": "center",
			"va": "bottom",
			"size": 14,
			"color": "black"
		}
		ax.annotate(**args)
	
		if row.Matched:
			xticklabels.append( (x, xlabel) )
	
	ax.set_xticks( list(map(lambda x:x[0], xticklabels)) )
	ax.set_xticklabels( list(map(lambda x:x[1], xticklabels)) )
	ax.ticklabel_format(axis="y", style="sci", scilimits=(6,6))
	ax.set_ylabel("# Reads")
	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/10)
	ax.legend(loc="upper left")
	
	sub_title = "\n({:,} matched/{:,} reads)".format(matched, total)
	ax.set_title( title + sub_title , size=18 )

	fig.tight_layout()
	return df, fig
	############################################################################
	
##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path, header=None, names=["Matched", "Mapped", "Reads"])
	df["Matched"] = df.Matched.replace({"MATCHED": True, "UNMATCHED": False})
	df["Mapped"] = df.Mapped.replace({"MAPPED": True, "UNMAPPED": False})
	
	ddf, plt = align_plot(df, "UP primer matching")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

