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
from scipy import optimize


from matplotlib.backends.backend_pdf  import PdfPages


def fitfunc (x, m, b):
    x=np.array(x)
    return m*x + b
pinit=[0.1, 0.1]

def linear_fit (x,y,y_err=None):
    pfinal, covar=optimize.curve_fit(fitfunc,x, y, p0=pinit, sigma=y_err,  maxfev=100000)
    return pfinal[0], np.sqrt(covar[0]), pfinal[1], np.sqrt(covar[1])

pp=PdfPages('histogram_sizes_substructures.pdf')

data=np.genfromtxt ("star_forming_sizes.dat")
size_vec_pix=data[:,0]
size_vec_kpc=data[:,1] # Actually in pc
nbins=25
cut=600 #pc
index=size_vec_kpc<cut
size_vec_kpc=size_vec_kpc[index]
size_vec_pix=size_vec_pix[index]


fig=plt.figure()
ax = fig.add_subplot(211)
n, bins, patches = ax.hist(size_vec_pix, nbins, normed=False, facecolor='green', alpha=0.75, label='')
plt.xlabel('Size (pixels)')
plt.ylabel('histogram')

ax = fig.add_subplot(212)
n, bins, patches = ax.hist(size_vec_kpc, nbins, normed=False, facecolor='red', alpha=0.75,label='')
plt.xlabel('Size (pc)')
plt.ylabel('histogram')
pp.savefig()

bincenters = 0.5*(bins[1:]+bins[:-1])

index=bincenters>200
x=bincenters[index]
n=n[index]

print "x, n: ", x, n
print "np.log(x),np.log(n): ", np.log(x),np.log(n)

m, m_err, b, b_err = linear_fit (np.log(x),np.log(n))
print "m, m_err: ", m, m_err
print "b, b_err: ", b, b_err


f=open('histogram_star_forming_sizes_pc.dat', 'w')
for i in range(len(n)):
    line="%g %g %g %g\n" %(bins[i],bincenters[i],bins[i+1],n[i])
    f.write(line)
f.close()


fig=plt.figure()
ax = fig.add_subplot(211)
n, bins, patches = ax.hist(size_vec_pix, nbins, normed=False, facecolor='green', alpha=0.75, label='')
plt.xlabel('Size (pixels)')
plt.ylabel('histogram')
plt.yscale('log',nonposy='clip')
plt.xscale('log')

ax = fig.add_subplot(212)
n, bins, patches = ax.hist(size_vec_kpc, bins=np.linspace(1,500,30), normed=False, facecolor='red', alpha=0.75,label='')
plt.xlabel('Size (pc)')
plt.ylabel('histogram')
plt.yscale('log',nonposy='clip')
plt.xscale('log')
plt.xlim([30,500])
pp.savefig()

bincenters = 0.5*(bins[1:]+bins[:-1])

print np.logspace(1, 500, 20)
print np.min(size_vec_kpc), np.max(size_vec_kpc)

for a,b in zip( bincenters, n):
    print a, b



pp.close()
