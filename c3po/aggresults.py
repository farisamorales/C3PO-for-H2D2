# Create the master result tables using the individual files from the fits.

import os
import re
import pandas as pd
import sed_config as sconf


def pull_direcs(direc):
    direcs = os.listdir(direc)
    indexes = []
    for i, d in enumerate(direcs):
        a = re.search("\w*\.\w+", d)
        if a: indexes.append(i)
    if len(indexes) == 0:
        return [direc + os.sep + direcs[i] for i in range(len(direcs))]
    for i in indexes[::-1]:
        direcs.pop(i)
    return [direc + os.sep + direcs[i] for i in range(len(direcs))]

def pull_fnames(direc):
    fnames = os.listdir(direc)
    keepers = []
    for i, fn in enumerate(fnames):
        a = re.search("(\w+\s\d+\.\w+)", fn)
        if a: keepers.append(fn)
    return keepers

def aggregate(direc):
    fnames = pull_fnames(direc)
    df = pd.read_csv(direc+os.sep+fnames[0])
    for fname in fnames[1:]:
        df2 = pd.read_csv(direc+os.sep+fname, index_col=False)
        df = pd.concat([df, df2])
    return df

direc = './Results/Params'
direcs = pull_direcs(direc)
for d in direcs:
    x = aggregate(d)
    x.to_csv(d+os.sep+'__MasterResults.csv', index=False)

direc = './Results/FluxRatios'
direcs = pull_direcs(direc)
for d in direcs:
    x = aggregate(d)
    x.to_csv(d+os.sep+'__MasterResults.csv', index=False)
