#!/usr/bin/python

import sys
from itertools import product
import numpy as np
from scipy.stats import entropy
import pandas as pd

hamming_path = sys.argv[1]
reads_path = sys.argv[2]
coords_path = sys.argv[3]
base_path = sys.argv[4]
max_distance = int(sys.argv[5])
max_matches = int(sys.argv[6])
max_entropy = float(sys.argv[7])

#hamming_path = "/camp/stp/babs/working/bahn/projects/tybulewiczv/micaela.sartoretti/PM21235/results/qc/25-06.ordered.hamming.csv"
#reads_path = "/camp/stp/babs/working/bahn/projects/tybulewiczv/micaela.sartoretti/PM21235/results/qc/25-06.reads_umis_per_barcode.csv"
#coords_path = "/camp/stp/babs/working/bahn/projects/tybulewiczv/micaela.sartoretti/PM21235/pucks/211025_06.csv"
#base_path = "tmp/matching"
#max_distance = 4
#max_matches = 5
#max_entropy = .5

#hamming_path = "results/qc/sample1.ordered.hamming.csv"
#hamming_path = "results/qc/sample1.shuffled.hamming.csv"
#reads_path = "results/qc/sample1.reads_umis_per_barcode.csv"
#hamming_path = "results/qc/sample2.shuffled.hamming.csv"
#hamming_path = "results/qc/sample2.ordered.hamming.csv"
#reads_path = "results/qc/sample2.reads_umis_per_barcode.csv"
#coords_path = "test/puck2.csv"
#base_path = "tmp/matching"
#max_distance = 4
#max_matches = 5
#max_entropy = .5

"""
The hamming distance data frame contains 4 colummns:
 * **SeqBarcode**: the illumina barcode
 * **Distance**: the smallest hamming distance that can be computed between an illumin barcode and one or several puck barcode
 * **Number**: the number of puck barcodes whose the hamming distance from the illumina barcode is the smallest that can be found
 * **PuckBarcode**: the list of puck barcodes (up to 10 items) whose the hamming distance from the illumina barcode is the smallest that can be found

The procedure includes 3 successive steps and is the following:
 1. If the distance is 0 then there is only one match, so we include it;
 2. If the distance is greater than 0 but there is only one match then we include it no matter what the distance is;
 3. If the match doesn't fall into the situations 1 or 2, then if the distance is lesser than or equal to {max_distance} and the number of matches is lesser than or equal to {max_matches}, then we take the match if the entropy of the reads in less than or equal than {max_entropy}.
"""

# minimum hamming distance for each, number of reads per barcodes and
# coordinates
hamming = pd.read_csv(hamming_path, header=None)
hamming.columns = ["SeqBarcode", "Distance", "Number", "PuckBarcode"]
hamming = hamming\
	.assign(PuckBarcode=lambda x: x.PuckBarcode.str.split(":"))\
	.explode("PuckBarcode")
reads = pd.read_csv(reads_path, header=None)
reads.columns = ["SeqBarcode", "Reads", "UMIs"]
coords = pd.read_csv(coords_path).rename(columns={"Barcode": "PuckBarcode"})

# default barcode
l = len(hamming.PuckBarcode.iloc[0])
bcd = "".join([ "X" for i in range(l) ])

# one big data frame
df = pd.merge(hamming, reads, how="outer")
df = pd.merge(df, coords, how="outer")
df = df\
	.assign(
		Sequenced=lambda x: np.where(x.SeqBarcode.isna(), "UNSEQUENCED", "SEQUENCED"),
		SeqBarcode=lambda x: np.where(x.SeqBarcode.isna(), bcd, x.SeqBarcode),
		Distance=lambda x: np.where(x.Distance.isna(), l, x.Distance),
		Number=lambda x: np.where(x.Number.isna(), 0, x.Number),
		Reads=lambda x: np.where(x.Reads.isna(), 0, x.Reads),
		UMIs=lambda x: np.where(x.UMIs.isna(), 0, x.UMIs)
	)

# working data frame
ddf = df.set_index("Sequenced").xs("SEQUENCED").reset_index(drop=True)

# perfect match
perfect = ddf.loc[ df.Distance == 0 ][ ["SeqBarcode", "PuckBarcode"] ]

####################################
# unperfect match but only one match

one = ddf
one = one.loc[ ~ one.SeqBarcode.isin(perfect.SeqBarcode) ]
one = one.loc[ ~ one.PuckBarcode.isin(perfect.PuckBarcode) ]
one = one.loc[ one.Number == 1 ][ ["SeqBarcode", "PuckBarcode"] ]

# hamming gives the number of matches for an illumina barcode, here we take
# the puck barcodes with only one match
indices = one\
	.groupby("PuckBarcode")\
	.apply(lambda x: x.shape[0] == 1)\
	.replace({False:np.nan})\
	.dropna()\
	.index

one = one.loc[ one.PuckBarcode.isin(indices) ]

##########
# the rest

################
def select(df):#
################
	values = df.Reads / df.Reads.sum()
	#values = df.UMIs / df.UMIs.sum()
	entrop = entropy(values, base=df.shape[0])

	if entrop <= max_entropy:
		status = "MATCHED"
		barcode = df.sort_values("Reads", ascending=False).SeqBarcode.iloc[0]
	else:
			status = "UNMATCHED"
			barcode = bcd
	
	return pd.Series([status, barcode], index=["Matched", "SeqBarcode"])
	############################################################################

more = ddf
more = more.loc[ ~ more.SeqBarcode.isin(perfect.SeqBarcode) ]
more = more.loc[ ~ more.PuckBarcode.isin(perfect.PuckBarcode) ]
more = more.loc[ ~ more.SeqBarcode.isin(one.SeqBarcode) ]
more = more.loc[ ~ more.PuckBarcode.isin(one.PuckBarcode) ]
more = more.loc[ more.Distance <= max_distance ]
more = more.loc[ more.Number <= max_matches ]

# mainly for the fake data, because on real data this always exists
if more.shape[0] > 0:
	more = more.groupby("PuckBarcode").apply(select)
	more = more.loc[ more.Matched == "MATCHED" ]
	more = more.reset_index().drop("Matched", axis=1)
	more = more.drop_duplicates("SeqBarcode", keep=False)
else:
	more = more[ ["SeqBarcode", "PuckBarcode"] ]

################
# final matching
matching = pd.concat([perfect, one, more])
matching["Matched"] = "MATCHED"

assert matching.SeqBarcode.unique().size == matching.shape[0], "Error illumina count is wrong"
assert matching.PuckBarcode.unique().size == matching.shape[0], "Error puck count is wrong"

d = pd.merge(df, matching, how="left")
d = d.assign(Matched=lambda x: np.where(x.Matched.isna(), "UNMATCHED", x.Matched))

#########
# metrics

n_seq = hamming.SeqBarcode.unique().size
n_puck = coords.PuckBarcode.unique().size
n_match = matching.shape[0]
n_perfect = perfect.shape[0]
n_one = one.shape[0]
n_more = more.shape[0]

metrics = pd.DataFrame.from_records([
	{"Metric": "Sequencing barcodes", "Value": n_seq},
	{"Metric": "Puck barcodes", "Value": n_puck},
	{"Metric": "Matches", "Value": n_match},
	{"Metric": "Distance 0", "Value": n_perfect},
	{"Metric": "Distance >0 and 1 match", "Value": n_one},
	{"Metric": f"Distance <{max_distance+1} and <{max_matches+1} match", "Value": n_more}
])

######
# save

final = d.loc[ d.Matched == "MATCHED" ]
final = final[ ["SeqBarcode", "PuckBarcode", "x", "y", "Distance", "Number", "Reads" ] ]
final.to_csv(f"{base_path}.csv", index=False)
d.to_csv(f"{base_path}.values.matching.csv", index=False)
metrics.to_csv(f"{base_path}.metrics.matching.csv", index=False)
matching.drop("Matched", axis=1).to_csv(f"{base_path}.map.matching.csv", index=False)

