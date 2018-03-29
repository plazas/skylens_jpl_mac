#!/usr/bin/env python

import numpy as np
import os
import sys
import math
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.font_manager as fm


import astropy.io.fits as pf

data=pf.open(sys.argv[1])[0].data
flux = []

pp=PdfPages (sys.argv[3])


N,M= data.shape
print data.shape


sum_data=np.sum(data)
print "sum_data", sum_data



# create data vector and unique index vector.
index_dict={}
k=0
index_vec=[]
for i in range(N):
    for j in range(M):
        #print "data: ", data[i,j]
        f=data[i,j]
        flux.append(f)
        index_dict[f]=[i,j,0]
        index_vec.append(k)
        k+=1
index_vec=np.array(index_vec)

### Histogram:

fig=plt.figure()
plt.hist(flux, bins=50, normed=True)
pp.savefig()

fig=plt.figure()
plt.hist(flux, bins=100, normed=True, cumulative=True)
pp.savefig()


sorted=np.argsort(flux)
flux=np.array(flux)
flux=flux[sorted]

index_vec=np.array(index_vec)
index_vec=index_vec[sorted]
cumu_prob=np.arange (len(flux))/ float(len(flux))

print "flux, cumu_prob ", flux, cumu_prob

fig=plt.figure()
plt.plot(flux,cumu_prob , 'm.')
pp.savefig()
pp.close()


for k,f in enumerate(flux):
    index_dict[f][2]=k


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

fluxes_sampled=[]
counter=0
n_samples=30
for u in np.random.uniform(0.95, 1.,n_samples):
    near_index, near_value = find_nearest (cumu_prob, u )
    print "random number, near index, near value: ", u, near_index, near_value
    fluxes_sampled.append(flux[near_index])

fluxes_sampled=np.array(fluxes_sampled)

final_data=np.zeros_like(data)


for f in fluxes_sampled:
    l=index_dict[f][0]
    m=index_dict[f][1]
    final_data[l,m]=data[l,m]+10000


hdu = pf.PrimaryHDU(final_data)
hdu.writeto (sys.argv[2], clobber=True)


sys.exit()






fig=plt.figure()
plt.plot( index_vec, flux, 'm.')
pp.savefig()

print "np.sum(flux)", np.sum(flux)

import esutil as eu
b=eu.stat.Binner(index_vec, flux)
#b.dohist(nperbin=10)
b.dohist(nbin=40)
b.calc_stats()
x=b['xmean']
y=b['ymean']
yerr=b['yerr']
print "x", x
print "y", y
print "yerr", yerr



fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot( x, y, 'm.')
pp.savefig()

delta_x=x[2]-x[1]
print "np.diff: ", np.diff(x)
print "np.sum(y)", np.sum(y*delta_x)


#cum=np.zeros_like (flux)
c=[0.]
y*=delta_x
for i in range(1,len(x)):
    print c[i-1] + y[i]
    c.append( (c[i-1] + y[i]))

print "cumulative: ", c

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( x, c)
pp.savefig()




pp.close()

pp.close()

sys.exit()

"""
cumulative=np.cumsum(flux)

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( index_vec, cumulative)
pp.savefig()



fig=plt.figure()
ax=fig.add_subplot(111)
plt.hist( flux, bins=len(flux), normed=1)
pp.savefig()

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( index_vec, flux, 'm.')
pp.savefig()


import esutil as eu
b=eu.stat.Binner(index_vec, flux)
#b.dohist(nperbin=10)
b.dohist(nbin=70)
b.calc_stats()
x=np.array(b['xmean'])
y=np.array(b['ymean'])
yerr=np.array(b['yerr'])

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( x, y, 'm.')
pp.savefig()

cum_values = np.zeros(x.shape)
cum_values[1:] = np.cumsum(y[1:]*np.diff(x))

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( x, cum_values, 'c.')
pp.savefig()
pp.close()

sys.exit()
"""

sorted=np.argsort(flux)
flux=np.array(flux)
flux=flux[sorted]
index_vec=index_vec[sorted]
print min(flux), max(flux)
print flux[0], flux[-1]

for f,k in zip(flux, index_vec):
    index_dict[f][2]=k

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( index_vec, flux, 'm.')
pp.savefig()

#cum=np.zeros_like (flux)
c=[flux[0]]

for i in range(1,len(flux)):
    c.append(c[i-1] + flux[i])

fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( index_vec, c, 'r.')
pp.savefig()


fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot( c, flux, 'b.')
pp.savefig()


fig=plt.figure()
ax=fig.add_subplot(111)
plt.plot(flux, c, 'b.')
pp.savefig()

#pp.close()


#pp.close()


import scipy.interpolate as interpolate
 
#def inverse_transform_sampling(data, n_bins=50, n_samples=10000):
#    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
#    print "bin_edges", bin_edges
#    cum_values = np.zeros(bin_edges.shape)
#    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))

#    fig=plt.figure()
#    ax=fig.add_subplot(111)
#    plt.plot(  bin_edges, cum_values, 'c.')
#    pp.savefig()
    
    
#    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
#    r = np.random.rand(n_samples)
#   return inv_cdf(r)


#positions = inverse_transform_sampling (flux, n_bins=len(flux)-1)

n_samples=40
#inv_cdf = interpolate.interp1d(cum, flux)
#r = np.random.rand(n_samples)
#sample=inv_cdf(r)


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

fluxes_sampled=[]
for u in np.random.rand(n_samples):
    near_index, near_value = find_nearest (c, u )
    print "u, near_value, index, cum[near_index], flux[near_index], index_dict[flux[near_index]]: ", u, near_value, near_index, c[near_index], flux[near_index], index_dict[flux[near_index]]
    fluxes_sampled.append(flux[near_index])

fluxes_sampled=np.array(fluxes_sampled)

#sample=inverse_transform_sampling(flux, n_bins=len(flux))
#print "sample", sample

fig=plt.figure()
ax=fig.add_subplot(111)
plt.hist( fluxes_sampled, color='r', bins=len(fluxes_sampled), normed=1)
#plt.hist( flux, color='b', bins=len(flux), normed=1, alpha=0.5)
pp.savefig()


#fig=plt.figure()
#ax=fig.add_subplot(111)
#plt.plot(flux, index, 'r.')
#pp.savefig()

#pp.close()
#sys.exit()


#index_func=interpolate.interp1d (flux, index)
#positions=index_func (sample)

#positions=positions.astype(int)
#print min(positions), max(positions), "positions"


final_data=np.zeros_like(data)


for f in fluxes_sampled:
    l=index_dict[f][0]
    m=index_dict[f][1]
    final_data[l,m]=data[l,m]+10000


hdu = pf.PrimaryHDU(final_data)
hdu.writeto (sys.argv[2], clobber=True)


pp.close()
