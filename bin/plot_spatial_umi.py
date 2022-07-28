#!/usr/bin/env python3

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scanpy import read_10x_mtx

import matplotlib as mpl
mpl.rc('font', size=16)

#####################################
def spatial_plot(df, column, title):#
#####################################

	fig = Figure(figsize=(8,8))
	ax = fig.add_subplot(111)
	c = ax.scatter(
		df.x, df.y,
		c=df[column],
		s=1,
		cmap="viridis_r",
		norm=mpl.colors.Normalize(0, np.percentile(df[column], 95), clip=True)
		)
	c.set_rasterized(True)
	ax.set_xlabel("X")
	ax.set_ylabel("Y")
	ax.axis("equal")
	ax.set_title(title)
	fig.colorbar(c, ax=ax)
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	mtx_path = sys.argv[1]
	coords_path = sys.argv[2]
	base_path = sys.argv[3]
	
	adata = read_10x_mtx(mtx_path)
	coords = pd.read_csv(coords_path).rename(columns={"PuckBarcode": "Barcode"})

	umis = adata\
		.to_df()\
		.apply(sum, axis=1)\
		.to_frame()\
		.reset_index()\
		.rename(columns={"index": "Barcode", 0: "UMIs"})

	df = pd.merge(umis, coords)

	plt = spatial_plot(df, "UMIs", "UMIs per bead")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

