#coding=utf-8
#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pyplot as plt
import esutil as eu
from scipy import ndimage
import sys
import matplotlib.cm as cm
import galsim
import subprocess as S
import astropy.io.fits as pf
import os
from matplotlib.colors import LogNorm


from matplotlib.backends.backend_pdf  import PdfPages


"""
Goes to NED website, get's redshifts and velocities for al the galaxies in Frei's catalog. 
Prints file 'z_distances_frei.dat': nXXXX | median(z) | dA (using Planck) kPc | r (v/H0) kPc. 
Then use this file in other code that finds star forming regions to convert their angular size to physical size.
"""

frei_images_path="/Users/amalagon/STRONG_LENSING/HST_SL_SIMS/data/frei_galaxy_catalog"
cmd="ls %s/*.fits" %frei_images_path
list=S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

c=eu.cosmology.Cosmo (H0=67.8, flat=True, omega_m=0.315)

def main (argv):
    out=open('z_distances_frei.dat','w')
    for file in list:
        #print "file: ", file
        root=file.split('/')[-1].split('_')[0][1:]
        print "root: ", root
        cmd="wget --no-check-certificate 'https://ned.ipac.caltech.edu/cgi-bin/datasearch?objname=NGC+%s&search_type=Redshifts&zv_breaker=30000.0&of=table' -O temp1.txt " %root
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        cmd="< temp1.txt awk '/z=/ {print $8}' > temp.dat"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        cmd="< temp1.txt awk '/z=/ {print $5}' > temp2.dat"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        
        #Angular diameter distance
        data= np.genfromtxt ("temp.dat")
        print "REDSHIFTS from NED: ", data
        z=np.median ( np.nan_to_num ( data ) )
        
        # hubble's Law for local universe
        data= np.genfromtxt ("temp2.dat")
        print data
        v=np.median ( np.nan_to_num ( data ) )
        r=(v/67.8)*1e3   # distance from Hubble's law in kPc. Should be similar to dA for low-z galaxies.
        da= (c.Da (0.0, z))*1e3  #kPc (Erin's routine returns Mpc)
        line="n%s %g %g %g \n" %(root, z, da, r)
        print "REDSHIFT, Angular diameter distance: ",  line
        out.write(line)
        cmd="rm temp1.txt temp.dat"
    out.close

if __name__ == "__main__":
    import pdb, traceback
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
