#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

####################
def error_plot(df):#
####################
	"""Plots histgram of ambiguous hits. The argument data frame has to contain
	at least 2 columns named "Base" and "Errors". The function returns a Figure
	object."""

	fig = Figure(figsize=(8,8))
	ax = fig.add_subplot(111)

	ax.bar(
		df.Base, df.Errors,
		width=1,
		color="skyblue",
		edgecolor="black",
	)

	n = df.Base.max()
	
	ax.set_xticks(range(1, n+1))
	xticklabels = []
	for i in range(1, n+1):
		if i % 2 == 1:
			xticklabels.append(str(i))
		else:
			xticklabels.append("")
	ax.set_xticklabels(xticklabels)
	
	ax.set_xlabel("Base position")
	ax.set_ylabel("Number of mismatches")
	ax.set_title(f"Histogram of the errors in matched barcodes")
	
	fig.tight_layout()

	return fig
	############################################################################

#######################
def error(seq1, seq2):#
#######################
	d = {}
	for i, (b1, b2) in enumerate( zip(seq1, seq2) , start=1 ):
		d[i] = 0 if b1 == b2 else 1
	return pd.Series(d)

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]
	
	df = pd.read_csv(csv_path)
	
	errors = df\
		.apply(lambda x: error(x.PuckBarcode, x.SeqBarcode), axis=1)\
		.apply(sum)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Base", 0: "Errors"})
	
	plt = error_plot(errors)
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

