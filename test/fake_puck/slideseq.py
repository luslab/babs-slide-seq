#!/usr/bin/python

import numpy as np
import pandas as pd
from skimage.io import imread
import matplotlib.pyplot as plt

X = imread("slideseq.tif")
zeros = np.where( X == 0 )

R = np.array(
    [np.cos(3*np.pi/2), -np.sin(3*np.pi/2), np.sin(3*np.pi/2), np.cos(3*np.pi/2)]
    ).reshape((2,2))

Y = np.matmul(R, np.array([zeros[0], zeros[1]]))
df = pd.DataFrame(Y.T)
df.columns = ["x", "y"]
df.to_csv("slideseq.csv", header=False, index=False, float_format="%.15f")

plt.scatter(df.x, df.y, s=.1)
plt.show()

