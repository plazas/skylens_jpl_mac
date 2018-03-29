#coding=utf-8
#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
import subprocess as S
import astropy.io.fits as pf
import os
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf  import PdfPages
from scipy import optimize

import matplotlib.font_manager as fm
import matplotlib.patches as patches


pp=PdfPages("F_distribution.pdf")
prop = fm.FontProperties(size=7)


Ftot=1.0 ## Always 1 for a PDF.
alpha=4.0
m=1-alpha
Fmin=10
Fmax=1e2

A=Ftot*m/(Fmax**m - Fmin**m)


u=np.random.uniform(0.,1.,10000)*Ftot
x = ((Fmax**m - Fmin**m)*u/Ftot + Fmin**m)**(1./m)


print len(u)
print len(x)
print min(x), max(x)



### Plot the power law

F=np.linspace(Fmin, Fmax, 1000)
p= A*F**-alpha

fig=plt.figure()
ax=fig.add_subplot(111)
ax.errorbar (F, p,yerr=None, fmt='r-', alpha=0.8, label='power law theory', markersize=3)
n, bins, patches = ax.hist(x, 50, normed=True, facecolor='green', alpha=0.75, label='F dist.')

#ax.errorbar (F, x,yerr=None, fmt='b-', alpha=0.8, label='power law', markersize=3)

plt.xlabel(r"F")
plt.ylabel('P')
plt.yscale('log',nonposy='clip')
plt.xscale('log')
plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
pp.savefig()


### Plot the cumulative

C=Ftot*(F**m - Fmin**m)*1./(Fmax**m - Fmin**m)
fig=plt.figure()
ax=fig.add_subplot(111)
ax.errorbar (F, C,yerr=None, fmt='g-', alpha=0.8, label='power law cumulative', markersize=3)
plt.xlabel(r"F")
plt.ylabel('P')
plt.yscale('log',nonposy='clip')
plt.xscale('log')
plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
pp.savefig()

pp.close()



