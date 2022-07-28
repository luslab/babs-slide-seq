#!/usr/bin/env python3
# coding: utf-8

import sys
import pandas as pd
from random import sample

path = sys.argv[1]
basepath = sys.argv[2]

df = pd.read_csv(path)
ordered = df.Barcode
shuffled = df.Barcode.apply(lambda x: "".join(sample(x, len(x))))
ordered.to_csv(f"{basepath}.ordered.txt", index=False, header=False)
shuffled.to_csv(f"{basepath}.shuffled.txt", index=False, header=False)
