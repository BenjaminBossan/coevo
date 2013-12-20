# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Plotting results from julia-generated csv's

# <headingcell level=2>

# Import packages

# <codecell>

import csv
import numpy as np
import matplotlib
import pylab as pl

# <headingcell level=3>

# dictionary containing the name of the strats

# <codecell>

stratDict = {"strat1": "indiv. learning",
             "strat2": "conformism",
             "strat3": "opp. indiv. learning",
             "strat4": "opp. conformism",
             "strat5": "In Doubt Conform",
             "strat6": "Imitate The Wealthiest",
             "strat7": "scoring PBSL, [4/-1]",
             "strat8": "scoring PBSL, [1/0]",
             "strat9": "averaging PBSL",
             "strat10": "avrg. PBSL w/ po-conf t-o"}

# <headingcell level=3>

# import data for plotting behavior of 1 generation

# <codecell>

fileName = "../julia/coevo/1gen01.csv"
d1 = []
i=0
with open(fileName, 'rb') as f:
    for row in csv.reader(f, delimiter = ";"):
        if i==0:
            header1 = row
        else:
            try:
                d1.append([float(d) for d in row])
            except:
                print i
                break
        i+=1
d1 = np.array(d1)

# <headingcell level=3>

# Plot

# <codecell>

pl.plot(d1[:, 0] - d1[:, 1] + .5, 'k-')
for col in np.transpose(d1[:, 2:]):
    pl.plot(col)
pl.ylim([0, 1])
legend1 = [stratDict[x] for x in header1[2:]]
legend1.insert(0, "pA - pB")
pl.legend(legend1)
pl.show()

# <headingcell level=3>

# import data for plotting evolution of strategy frequencies

# <codecell>

fileName2 = "../julia/coevo/evo01.csv"
d2 = []
i=0
with open(fileName2, 'rb') as f:
    for row in csv.reader(f, delimiter = ";"):
        if i==0:
            header2 = row
        else:
            try:
                d2.append([float(d) for d in row])
            except:
                print i
                break
        i+=1

# <codecell>

d2=np.array(d2)

# <headingcell level=3>

# Plot

# <codecell>

for col in np.transpose(d2):
    pl.plot(col)
pl.ylim([0, 1])
pl.legend([stratDict[x] for x in header2])
pl.show()

# <codecell>


