#coding=utf-8
#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import esutil
from scipy import ndimage
import sys
import matplotlib.cm as cm
import galsim
import subprocess as S
import astropy.io.fits as pf
import os
from matplotlib.colors import LogNorm

from scipy import stats
#from scipy import stsci
from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.font_manager as fm


pp=PdfPages('stamps.pdf')

file1="./obj_22245_no_knots_F475W.fits"
file2="./obj_22245_with_knots_F475W.fits"

data1=pf.open(file1)[0].data
data2=pf.open(file2)[0].data

factor1=0.01
factor2=2.5


mean_data1 = np.mean (data1)
mean_data2 = np.mean (data2)

print "mean data1: ", mean_data1
print "mean data2: ", mean_data2

fig=plt.figure()
ax=fig.add_subplot (211)
ax.imshow(data1, cmap=cm.gray, origin="lower", interpolation='nearest', vmin=mean_data1 - mean_data1/factor1, vmax=mean_data1 + mean_data1/factor1)
ax=fig.add_subplot (212)
ax.imshow(data2, cmap=cm.gray, origin="lower", interpolation='nearest', vmin=mean_data2 - mean_data2/factor2, vmax=mean_data2 + mean_data2/factor2)

plt.tight_layout()
pp.savefig()
pp.close()
