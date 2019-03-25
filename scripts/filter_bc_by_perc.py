#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

infile = sys.argv[1]
perc_spec = int(sys.argv[2])

df = pd.read_table(infile, sep="\t", header=None, names = ['count','bc'])
df[['count']] = df[['count']].astype(int)
df = df[df['count'] <= np.percentile(df['count'],perc_spec)]

df = df[['bc']]
df.to_csv("bcs_include_perc" + str(perc_spec) + ".txt", header=False, index=False)

