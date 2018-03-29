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

from astropy.convolution import convolve, Box2DKernel, Gaussian2DKernel
import esutil as eu

from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from scipy import optimize



def fitfunc (x, m, b):
    x=np.array(x)
    return m*x + b
pinit=[0.1, 0.1]

def linear_fit (x,y,y_err=None):
    pfinal, covar=optimize.curve_fit(fitfunc,x, y, p0=pinit, sigma=y_err,  maxfev=100000)
    return pfinal[0], np.sqrt(covar[0]), pfinal[1], np.sqrt(covar[1])



def draw_sersic_model (mag=-10, r=1, pixel_scale=1, q=0.2, beta=0.0, n=1, xc=0, yc=0, x_size=64, y_size=64):
    k=64
    big_fft_params = galsim.GSParams(maximum_fft_size=int(512*k))

    flux_gal=10**(-0.4*mag)  # Lensed documentation: L_tot = 10**(-0.4m)
    size_gal=r         # In arcseconds
    q_gal=q
    beta_gal= beta
    n=n
    #random_seed=1534225
    #rng = galsim.BaseDeviate(random_seed)
    xc=xc
    yc=yc
    
    gal=galsim.Sersic(n=n, half_light_radius=size_gal, flux=flux_gal, gsparams=big_fft_params).shear(q=q_gal, beta=beta_gal*galsim.degrees )
    base_size_x=x_size
    base_size_y=y_size
    im=gal.drawImage(image=galsim.ImageF(base_size_x, base_size_y),scale=pixel_scale, offset=(xc, yc))
    return im.array


def draw_bulge_plus_disk_model (mag_b=-10, r_b=1, q_b=0.2, beta_b=0.0, n_b=1, mag_d=-10, r_d=1, q_d=0.2, beta_d=0.0, n_d=1, pixel_scale=1, xc=0, yc=0, x_size=64, y_size=64):
    k=64
    big_fft_params = galsim.GSParams(maximum_fft_size=int(512*k))
    
    flux_b=10**(-0.4*mag_b)
    bulge = galsim.Sersic(n=n_b, half_light_radius=r_b, flux=flux_b, gsparams=big_fft_params).shear(q=q_b, beta=beta_b*galsim.degrees)

    flux_d=10**(-0.4*mag_d)
    disk = galsim.Sersic(n=n_d, half_light_radius=r_d, flux=flux_d, gsparams=big_fft_params).shear(q=q_d, beta=beta_d*galsim.degrees)
    
    gal=disk+bulge
    
    xc=xc
    yc=yc
    
    base_size_x=x_size
    base_size_y=y_size
    im=gal.drawImage(image=galsim.ImageF(base_size_x, base_size_y),scale=pixel_scale, offset=(xc, yc))
    return im.array




def get_number (cmd):
    return float(S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()[0])

def histeq(im,nbr_bins=256):
    
    #get image histogram
    imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
    cdf = imhist.cumsum() #cumulative distribution function
    cdf = 255 * cdf / cdf[-1] #normalize
    #use linear interpolation of cdf to find new pixel values
    im2 = np.interp(im.flatten(),bins[:-1],cdf)
    return im2.reshape(im.shape), cdf


to_rad=np.pi*1./180

def is_in_ellipse (x,y, xc, yc, a, b, theta):
    if ( ((x-xc)*np.cos(theta*to_rad) + (y-yc)*np.sin(theta*to_rad))**2/a**2+ ((x-xc)*np.sin(theta*to_rad)-(y-yc)*np.cos(theta*to_rad))**2/b**2) <= 1 :
        return 0.0
    else:
        return 1.0




#### PLOTS
#### Do the plotting here
plt.minorticks_on()
#plt.tight_layout()

### We do not have matplotlib 1.1, with the 'style' package. Modify the matplotlibrc file parameters instead
import matplotlib as mpl
mpl.rc('lines', linewidth=1, color='black', linestyle='-')
mpl.rc('font', family='serif',weight='normal', size=10.0 )
mpl.rc('text',  color='black', usetex=False)
mpl.rc('axes',  edgecolor='black', linewidth=1, grid=False, titlesize=14, labelsize=14, labelweight='normal',labelcolor='black')
mpl.rc('axes.formatter', limits=[-4,4])
mpl.rcParams['xtick.major.size']=14
mpl.rcParams['xtick.minor.size']=8
mpl.rcParams['xtick.major.pad']=8
mpl.rcParams['xtick.minor.pad']=8
mpl.rcParams['xtick.labelsize']= '12'
mpl.rcParams['xtick.minor.width']= 1.0
mpl.rcParams['xtick.major.width']= 1.0
mpl.rcParams['ytick.major.size']=14
mpl.rcParams['ytick.minor.size']=8
mpl.rcParams['ytick.major.pad']=8
mpl.rcParams['ytick.minor.pad']=8
mpl.rcParams['ytick.labelsize']= '12'
mpl.rcParams['ytick.minor.width']= 1.0
mpl.rcParams['ytick.major.width']= 1.0
mpl.rc ('legend', numpoints=1, fontsize='14', shadow=False, frameon=False)






print "Starting: "
output_directory_path="/Users/amalagon/STRONG_LENSING/SKYLENS3_SIMULATIONS/FREI_GUNN_GALAXIES_LENSED/OUTPUT_LENSED_BULGE_PLUS_DISK/"
cmd="ls -d %s*/" %output_directory_path
list_output_dirs=S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
frei_images_path="/Users/amalagon/STRONG_LENSING/SKYLENS3_SIMULATIONS/data/frei_galaxy_catalog/"
frei_distances_path="/Users/amalagon/STRONG_LENSING/SKYLENS3_SIMULATIONS/FREI_GUNN_GALAXIES_LENSED/z_distances_frei.dat"  # n2403 0.000437 (z) 1931.25 (dA, kpc) 1843.66 (hubble d, kpc)

arcsec_to_rad = 1./206264.81
prop = fm.FontProperties(size=7)
sigma_cut=10
#Read properties of galaxies
prop=np.genfromtxt("./properties_of_galaxies_name_type.txt") #NGC   Obs.     T  alpha        delta        B_T    v      D    I   p.a.
NGC=prop[:,0]
T=prop[:,1]
#print NGC
#print T

root_dir_frei_stamps="/Users/amalagon/STRONG_LENSING/SKYLENS3_SIMULATIONS/FREI_GUNN_GALAXIES_LENSED/OUTPUT_LENSED_BULGE_PLUS_DISK/"
lista=['n2403_pr', 'n3938_lr', 'n2715_lr', 'n2903_pr', 'n3319_pr', 'n3344_lr','n3368_lr',\
       'n3486_lr', 'n3631_lr', 'n3810_lr', 'n3893_lr', 'n3953_lr',\
       'n4030_lr', 'n4136_lr', 'n4178_pr', 'n4254_pr', 'n4258_pr', 'n4303_pr', 'n4321_pr', \
'n4501_pr', 'n4527_pr', 'n4535_pr', 'n4654_pr', 'n5055_pr', 'n5248_lr', 'n5364_lr', 'n5371_lr', \
'n5669_lr', 'n5701_lr', 'n6015_lr', 'n6118_lr', 'n6384_lr']



gauss_filter_size={'n2403_pr':3.5, 'n3938_lr':4.0, 'n2715_lr': 3, 'n2903_pr': 3.5, 'n3319_pr': 3.5, \
'n3344_lr': 3.5, 'n3368_lr': 4.0, 'n3486_lr': 10.0, 'n3596_lr': 8.0, 'n3631_lr': 15.0, 'n3672_lr':10.0, \
'n3810_lr':10.0, 'n3893_lr': 5.0, 'n3953_lr': 5.0, 'n4030_lr': 8.0, 'n4136_lr': 4.0, 'n4178_pr': 5.0, \
'n4192_pr': 7.0, 'n4254_pr': 10.0, 'n4258_pr': 10.0, 'n4303_pr': 10.0, 'n4321_pr':11.0, 'n4501_pr': 10.0, \
    'n4527_pr': 10.0, 'n4535_pr': 10.0, 'n4654_pr': 8.0, 'n5055_pr': 10.0, 'n5248_lr': 10.0, 'n5364_lr': 10.0,\
'n5371_lr': 10.0, 'n5669_lr': 10.0, 'n5701_lr': 10.0, 'n6015_lr': 10.0, 'n6118_lr': 10.0, 'n6384_lr': 10.0}

factor_ellipse={'n2403_pr':0.45, 'n3938_lr':0.6 , 'n2715_lr':0.7, 'n2903_pr': 0.8, 'n3319_pr': 0.7, \
'n3344_lr': 0.95, 'n3368_lr':0.85, 'n3486_lr':0.8, 'n3596_lr': 0.4, 'n3631_lr': 0.5, 'n3672_lr': 1.0,\
'n3810_lr':0.55, 'n3893_lr': 1.1, 'n3953_lr': 1.0, 'n4030_lr': 0.6, 'n4136_lr': 0.6, 'n4178_pr': 0.6,\
'n4192_pr':0.8, 'n4254_pr': 0.8, 'n4258_pr': 0.75, 'n4303_pr': 0.3, 'n4321_pr':0.4, 'n4501_pr':0.9,
'n4527_pr': 0.9, 'n4535_pr': 0.35, 'n4654_pr': 0.65, 'n5055_pr': 0.7, 'n5248_lr': 0.5, 'n5364_lr': 0.7,\
'n5371_lr':0.6, 'n5669_lr': 0.9, 'n5701_lr': 0.9, 'n6015_lr': 1.0, 'n6118_lr': 0.9, 'n6384_lr': 0.6}

factor_contrast={'n2403_pr':1.0, 'n3938_lr':1.0 , 'n2715_lr':1.0, 'n2903_pr': 1.0, 'n3319_pr': 1.0, \
'n3344_lr': 1.4, 'n3368_lr':2.0,'n3486_lr':3.4, 'n3596_lr': 0.4, 'n3631_lr': 1.0, 'n3672_lr':0.25,
'n3810_lr':0.75, 'n3893_lr':0.5, 'n3953_lr': 1.0, 'n4030_lr':0.5, 'n4136_lr':0.75, 'n4178_pr': 0.8,\
'n4192_pr':0.7, 'n4254_pr': 0.5, 'n4258_pr':0.3, 'n4303_pr':0.3, 'n4321_pr':0.5, 'n4501_pr':0.5, \
'n4527_pr': 0.7, 'n4535_pr': 0.75, 'n4654_pr':0.8, 'n5055_pr': 0.9, 'n5248_lr': 0.4, 'n5364_lr': 1.5,\
'n5371_lr':0.8, 'n5669_lr': 1.2, 'n5701_lr': 1.7, 'n6015_lr': 0.5, 'n6118_lr': 1.5, 'n6384_lr': 1.5}

print "Number of stamps selected: ", len(lista)

type_dict={}
for a,b in zip(NGC,T):
    type_dict[a]=b

def main (argv):
    pp=PdfPages("out.pdf")
    f_aper_vec=[]
    f_auto_vec=[]
    size_vec_pix=[]
    size_vec_kpc=[]
    e_vec=[]
    chi2_vec=[]
    redshift_vec=[]
    dA_dict={}
    data=np.genfromtxt(frei_distances_path, dtype=[('gal','|S5'),('z','f8'), ('da','f8'), ('r_hubble','f8')])
    for i, gal in enumerate(data['gal']):
        dA_dict[gal]=np.abs(data['da'][i])

    for stamp in lista:
    #for stamp in ['n5364_lr']: #,'n2403_pr', 'n3938_lr']:
        dir=root_dir_frei_stamps+stamp+'/'
        print "dir: ", dir
        print type(dir)
        root=dir.split('/')[-2]
        print "root: ", root
        ngc=float(root.split('_')[0][1:])
        hd=pf.open(frei_images_path+root+".fits")[0].header
        #Select only spirals T in [0,9]
        T=int(type_dict[ngc])
        print "Morphological type T: ", T
        if not T in range(0,10):
            print "Skipping: T is not in [0,9] (i.e., not a spiral galaxy). "
            continue
        tel=hd['TELESCOP']
        print "tel: ", tel
        if tel == 'Palomar 1.5':
            pixel_scale = 1.19  #arcsec per pixel
        elif tel == 'Lowell 1.1':
            pixel_scale = 1.35
        else:
            print "Warning: Telescope is not Palomar nor Lewis. Setting plate scale to 1."
            pixel_scale=1.0
	    #read parameters of Sérsic model measured by LENSED
        file=dir+root+"-stdot.dat"
        cmd="< %s awk '/chi/ {print $4}'" %file
        if len(S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()) == 0:
            continue
        chi2=S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()[0]
        chi2=chi2[4:]  # some weird characters at the beginning of the string
        #chi2_vec.append(float(chi2))
        #if float(chi2) > 4:
        #    print "Skipping: chi2 larger than 4"
        #    continue


        print "Chi2: ", chi2
        #bulge
        cmd="< %s awk '/bulge.x/ {print $2}'" %file
        xc_b=get_number(cmd)
        cmd="< %s awk '/bulge.y/ {print $2}'" %file
        yc_b=get_number(cmd)
        cmd="< %s awk '/bulge.r/ {print $2}'" %file
        r_b=get_number(cmd)
        cmd="< %s awk '/bulge.mag/ {print $2}'" %file
        mag_b=get_number(cmd)
        cmd="< %s awk '/bulge.n/ {print $2}'" %file
        n_b=get_number(cmd)
        cmd="< %s awk '/bulge.q/ {print $2}'" %file
        q_b=get_number(cmd)
        cmd="< %s awk '/bulge.pa/ {print $2}'" %file
        beta_b=get_number(cmd)
        print "Bulge model parameters: ", xc_b, yc_b, r_b, mag_b, n_b, q_b, beta_b


        #disk
        cmd="< %s awk '/disk.x/ {print $2}'" %file
        xc_d=get_number(cmd)
        cmd="< %s awk '/disk.y/ {print $2}'" %file
        yc_d=get_number(cmd)
        cmd="< %s awk '/disk.r/ {print $2}'" %file
        r_d=get_number(cmd)
        cmd="< %s awk '/disk.mag/ {print $2}'" %file
        mag_d=get_number(cmd)
        cmd="< %s awk '/disk.n/ {print $2}'" %file
        n_d=get_number(cmd)
        cmd="< %s awk '/disk.q/ {print $2}'" %file
        q_d=get_number(cmd)
        cmd="< %s awk '/disk.pa/ {print $2}'" %file
        beta_d=get_number(cmd)
        print "Disk model parameters: ", xc_d, yc_d, r_d, mag_d, n_d, q_d, beta_d


        # Get model image

        xc = 0.5*(xc_b + xc_d)
        yc = 0.5*(yc_b + yc_d)
        
        print "Buldge centroid: ", xc_b, yc_b
        print "Disk centroid: ", xc_d, yc_d
        print "Mean centroid (from bulge + disk models): ", xc, yc

        tel=root.split('_')[1][0]
        if tel == 'p':
            pixel_scale=1.19
        else:
            pixel_scale=1.34
        x_size=hd['NAXIS1']
        y_size=hd['NAXIS2']
        x_offset = xc - x_size*1./2
        y_offset = yc - y_size*1./2
        #print "offsets: ", x_offset, y_offset


        #model_im = draw_sersic_model (mag=mag, r=r, pixel_scale=pixel_scale, q=q, beta=beta, n=n, xc=x_offset, yc=y_offset, x_size=x_size, y_size=y_size)


        model_im = draw_bulge_plus_disk_model (mag_b=mag_b, r_b=r_b, q_b=q_b, beta_b=beta_b, n_b=n_b, mag_d=mag_d, r_d=r_d, q_d=q_d, beta_d=beta_d, n_d=n_d, pixel_scale=pixel_scale, xc=x_offset, yc=y_offset, x_size=x_size, y_size=y_size)


        # Real data image
        print "Astronomical image: ", frei_images_path+root+".fits"
        data=pf.open(frei_images_path+root+".fits")[0].data
        
        ### MAKE A BAD PIXEL MASK OF CENTRAL REGION OF GALAXY
        bad=np.ones(data.shape)
        
        
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val=is_in_ellipse (i,j,xc_b,yc_b, r_b*factor_ellipse[stamp], q_b*(r_b*factor_ellipse[stamp]), beta_b)
                bad[j][i]=val
    
        #Write the bad pixel mask file to a fits file
        cmd="rm ./weight.fits"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        pf.writeto("./weight.fits", bad, output_verify='ignore')

        
        # Difference of original data and smoothed data
        data_smooth=convolve(data, Gaussian2DKernel( gauss_filter_size[stamp] ))
        diff=(data - data_smooth)

        #Write the diff file to a fits file, and then run Source Extractor, with the weight file from above
        cmd="rm diff.fits"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        pf.writeto("diff.fits", diff, output_verify='ignore')
        
        output_catalog="OUTPUT_SEXTRACTOR_BULGE_PLUS_DISK/%s_out.cat"%root
        cmd="mkdir -v OUTPUT_SEXTRACTOR_BULGE_PLUS_DISK"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()


        cmd="rm check.fits"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()


        sextractor_config="./sextractor_config_files/daofind_sex_%s.config" %stamp
        cmd="sex diff.fits -c %s -CATALOG_NAME %s" %(sextractor_config, output_catalog)
        print "Running SEXTRACTOR: "
        print cmd
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

        cmd="cp check.fits OUTPUT_SEXTRACTOR_BULGE_PLUS_DISK/%s_check.fits"%root
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        ## Open sextractor catalog and get the star forming regions
        cat=np.genfromtxt(output_catalog, dtype=[ ('FLUX_APER', 'f8'), ('FLUXERR_APER', 'f8'),('FLUX_AUTO', 'f8'),
                                                  ('FLUXERR_AUTO', 'f8'), ('FLUX_MAX', 'f8'), ('XWIN_IMAGE', 'f8'), ('YWIN_IMAGE', 'f8'),
                                                   ('ELLIPTICITY', 'f8'), ('FLUX_RADIUS ', 'f8'), ('CLASS_STAR','f8')] )
        print " len(np.atleast_1d(cat)): ",  len(np.atleast_1d(cat))
        if  len(np.atleast_1d(cat))==1: continue
        #Contents of the sextractor catalogs
        #   1 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
        #   2 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
        #   3 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
        #   4 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
        #   5 FLUX_MAX               Peak flux above background                                 [count]
        #   6 XWIN_IMAGE             Windowed position estimate along x                         [pixel]
        #   7 YWIN_IMAGE             Windowed position estimate along y                         [pixel]
        #   8 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
        #   9 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
        #  10 CLASS_STAR             S/G classifier output
        
        cat.sort(order='FLUX_RADIUS')
        try:
            x=cat['XWIN_IMAGE']
            y=cat['YWIN_IMAGE']
            r=cat['FLUX_RADIUS']
            f_aper=cat['FLUX_APER']
            f_auto=cat['FLUX_AUTO']
            e=cat['ELLIPTICITY']
        except IndexError:
            continue

        #distance_from_centroid = np.sqrt ( (x-xc)**2 + (y-yc)**2 )
        #index=(r>0.001) & (r < 20.0) & (e < 0.9) #& (distance_from_centroid > 5) # & (f_auto < 2e3) & (e < 0.5) # in pixel

        index=(x> 5) & (np.abs(y-data.shape[0]) > 10) & (np.abs(x-data.shape[1]) > 10) & (f_auto < 28000)
        r=r[index]
        f_aper=f_aper[index]
        f_auto=f_auto[index]
        e=e[index]
        x=x[index]
        y=y[index]


        print "len(r)", len(r)
        if len (r) == 0:
            print "len(r) == 0. Skipping. "
            continue
        #p=np.percentile (r, 40)  # choose the smaller objects found by Sextractor
        #print "Pixel cut in object radius: ", p
        total_flux_objects=f_aper.sum()
        #index=r<p #selected objects
        #print "len(index): ", len(index)
        #r=r[index]  #radius in pixels
        #f_aper=f_aper[index]/ total_flux_objects
        #f_auto=f_auto[index]
        #e=e[index]
        #x=x[index]
        #y=y[index]

        print "x: ", x
        print "y: ", y
        print "r: ", r
        print "f auto: ", f_auto
        size_vec_pix.extend(r)
        f_aper_vec.extend(f_aper)
        f_auto_vec.extend(f_auto)
        chi2_vec.append(float(chi2))
        e_vec.extend(e)


        dA = dA_dict[root[0:5]]
        r_kpc = dA*(r*pixel_scale*arcsec_to_rad)  #dA in kpc
        r_pc=r_kpc*1e3
        size_vec_kpc.extend(r_pc)
        print "Size in pixels; size in pc: ", r, r_pc


        #load the segmentation file
        ap=pf.open("check.fits")[0].data

        #fig=plt.figure()
        #ax=fig.add_subplot(221)
        #n, bins, patches = ax.hist(data.flatten(), 50, normed=False, facecolor='green', alpha=0.75, label='')
        #ax=fig.add_subplot(222)
        #n, bins, patches = ax.hist(model_im.flatten(), 50, normed=False, facecolor='green', alpha=0.75, label='')
        #ax=fig.add_subplot(223)
        #n, bins, patches = ax.hist(diff.flatten(), 50, normed=False, facecolor='green', alpha=0.75, label='')
        #ax=fig.add_subplot(224)
        #n, bins, patches = ax.hist(ap.flatten(), 50, normed=False, facecolor='green', alpha=0.75, label='')
        #pp.savefig()


        factor=factor_contrast[stamp]
        #data=ndimage.filters.gaussian_laplace(data,3)
        #data= convolve(data, Box2DKernel(3))
        mean_data = np.mean (data)
        print "mean data: ", mean_data
        fig=plt.figure()

        #ax=fig.add_subplot (121)  ### THIS IS FOR PAPER
        ax=fig.add_subplot (221)

        #data, d=histeq(data)

        ax.imshow(data, cmap=cm.gray, origin="lower", interpolation='nearest', vmin=mean_data - mean_data/factor, vmax=mean_data + mean_data/factor)
        #ax=plt.gca() #get the current axes
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        #circ = Circle((xc_b,yc_b), r_b*0.5, fill=False, color='green')
        #ax.add_patch(circ)


        #plt.colorbar()
        ax.set_title ("data", size=14)
        ax.tick_params(labelsize=12)

        #"""# FOR PAPER COMMENT THIS
        mean_model = np.mean (data_smooth)
        #print "data_smooth: ", data_smooth
        ax=fig.add_subplot (223)
        #ap, d=histeq(ap)
        #data_smooth=ap
        ax.imshow((data_smooth), origin="lower", cmap=cm.gray, interpolation='nearest', vmin=mean_model - mean_model/factor, vmax=mean_model + mean_model/factor)
        #ax.imshow((bad), origin="lower", cmap=cm.gray, interpolation='nearest')#, vmin=0, vmax=1)
        #ax.set_title ("model (chi2: %s)" %chi2, size=9)
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        #plt.colorbar()
        el=Ellipse( (xc_b,yc_b), 2*(r_b*factor_ellipse[stamp]), 2*(r_b*factor_ellipse[stamp])*q_b, angle=beta_b, fill=False, color='green')
        ax.add_patch(el)
        ax.set_title ("data smoothed", size=14)
        ax.tick_params(labelsize=2)
        #"""  # FOR PAPER COMMENT THIS

        mean_diff = np.mean (data)
        print "diff: ", mean_diff
        ax=fig.add_subplot (122)
        #diff,d=histeq(diff)
        ax.imshow((data), cmap=cm.gray, origin="lower", interpolation='nearest', vmin=mean_diff - mean_diff/factor, vmax=mean_diff + mean_diff/factor)
        #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
        #plt.colorbar(PCM, ax=ax)
        #plt.colorbar()
        for xx,yy in zip(x,y):
            circ = Circle((xx,yy),4, fill=False, color='red')
            ax.add_patch(circ)
        el=Ellipse( (xc_b,yc_b), 2*(r_b*factor_ellipse[stamp]), 2*(r_b*factor_ellipse[stamp])*q_b, angle=beta_b, fill=False, color='green')
        ax.add_patch(el)

        #ax.set_title ("diff", size=9)
        ax.set_title ("substructures", size=14)
        ax.tick_params(labelsize=12)



        """
        #mean_ap, scatter_ap, indices = eu.stat.sigma_clip (ap.flatten(), niter=10, nsig=sigma_cut, get_indices=True)
        mean_ap, scatter_ap= np.mean(ap), np.std(ap)
        ax=fig.add_subplot (224)
        print "max and min of diff: ", np.max(diff), np.min(diff)
        #ap,d=histeq(ap)
        ax.imshow(ap, cmap=cm.gray , origin="lower", interpolation='nearest' , vmin=-10, vmax=10)
        #for xx,yy in zip(x,y):
        #    circ = Circle((xx,yy),6, fill=False)
        #    ax.add_patch(circ)
        ax.set_title ("data - smoothed image", size=9)
        """
        plt.suptitle("%s"%("NGC"+root[1:5]), size=14)
        plt.tight_layout()

        pp.savefig()

    print " "
    print "Loop finished. "
    print " "

    size_vec_pix=np.array(size_vec_pix)
    size_vec_kpc=np.array(size_vec_kpc)
    f_aper_vec=np.array(f_aper_vec)
    f_auto_vec=np.array(f_auto_vec)
    chi2_vec=np.array(chi2_vec)
    e_vec=np.array(e_vec)
    
    print len(  size_vec_pix), len(  size_vec_kpc), len(f_aper_vec), len(f_auto_vec), len(chi2_vec)
    #index= (size_vec_kpc < 5000)
    print len(index)
    
    
    #size_vec_kpc=size_vec_kpc[index]
    #size_vec_pix=size_vec_pix[index]
    #f_aper_vec=f_aper_vec[index]
    #f_auto_vec=f_auto_vec[index]
    #chi2_vec=chi2_vec[index]
    
    
    #### Do sigma clipping in size
    import esutil as eu
    sigma_cut=3
    mean_size, std_size, index = eu.stat.sigma_clip (size_vec_kpc, niter=10, nsig=sigma_cut, get_indices=True, verbose=True)
    size_vec_kpc=size_vec_kpc[index]
    size_vec_pix=size_vec_pix[index]
    f_aper_vec=f_aper_vec[index]
    f_auto_vec=f_auto_vec[index]
    #chi2_vec=chi2_vec[index]
    
    print "size_vec_pix", len(size_vec_pix)
    print "size_vec_pc", len(size_vec_kpc)
    print "f_aper_vec", len(f_aper_vec)
    print "f_auto_vec", len(f_auto_vec)
    print "e_vec", len(e_vec)
    print "chi2_vec", len(chi2_vec)
    
    size_file=open("star_forming_sizes.dat", 'w')
    for p, x, flujo in zip(size_vec_pix, size_vec_kpc, f_auto_vec):
        line="%g %g %g \n"%(p,x,flujo)   # Units: pixels, parsecs! (not kpc)
        size_file.write(line)
    size_file.close()


    ntotal=len(size_vec_kpc)
    #mean_size=np.mean(size_vec_kpc)
    #std_size=np.std(size_vec_kpc)

    print "ntotal, mean_size, std_size  ", ntotal, mean_size, std_size

    fig=plt.figure()
    print "hola 1"
    ax = fig.add_subplot(211)
    n, bins, patches = ax.hist(size_vec_pix, 40, normed=False, facecolor='green', alpha=0.75)
    plt.xlabel('Size (pixels)')
    plt.ylabel('Frequcny')

    ax = fig.add_subplot(212)
    print "hola 2"
    n, bins, patches = ax.hist(size_vec_kpc, 40, normed=False, facecolor='red', alpha=0.75,label="Total number (after sigma-%g clipping): %g \n Mean: %g \n Std: %g" %(sigma_cut, ntotal,mean_size, std_size))
    plt.xlabel('Size (pc)')
    plt.ylabel('Frequency')
    #plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.suptitle("Total number: %g. Mean: %g pc. Std: %g pc." %(ntotal,mean_size, std_size))
    plt.tight_layout()
    pp.savefig()


    fig=plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = ax.hist(size_vec_kpc, 25, normed=False, facecolor='red', alpha=0.75,label="Total number: %g \n Mean: %g \n Std: %g" %(ntotal,mean_size, std_size))
    plt.xlabel('Half-light radius (pc)', size=16)
    plt.ylabel('Frequency', size=16)
    #plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)
    plt.tight_layout()
    pp.savefig()

    ##### Fit to a power law 
    bincenters_size = 0.5*(bins[1:]+bins[:-1])
    cut_size = 200
    index=bincenters_size>cut_size
    x_size=bincenters_size[index]
    n_size=n[index]

    print "x_size, n_size: ", x_size, n_size
    print "np.log(x_size),np.log(n_size): ", np.log(x_size),np.log(n_size)

    m_size, m_err_size, b_size, b_err_size = linear_fit (np.log(x_size),np.log(n_size))
    print "m_size, m_err_size: ", m_size, m_err_size
    print "b_size, b_err_size: ", b_size, b_err_size
    print " "

   
    
    mean_f_auto=np.mean(f_auto_vec)
    std_f_auto=np.std(f_auto_vec)

    print " mean_f_auto, std_f_auto ",  mean_f_auto, std_f_auto
    fig=plt.figure()
    #ax = fig.add_subplot(211)
    #print "hola 3"
    #n, bins, patches = ax.hist(chi2_vec, 50, normed=False, facecolor='blue', alpha=0.75, label='')
    #plt.xlabel('Reduced chi2')
    #plt.ylabel('histogram')
    print "hola 3"

    ax = fig.add_subplot(211)
    n, bins, patches = ax.hist(f_auto_vec, 15, normed=False, facecolor='yellow', alpha=0.75, label="Mean: %g \n Std: %g" %(mean_f_auto, std_f_auto))
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('F_AUTO')
    plt.ylabel('Frequency')
    #plt.legend(loc='upper left', fancybox=True, ncol=1, numpoints=1, prop = prop)


    print "Calculating power law parameters for luminosity: "

    bincenters = 0.5*(bins[1:]+bins[:-1])

    index=bincenters>2000
    x=bincenters[index]
    n=n[index]

    print "x, n: ", x, n
    print "np.log(x),np.log(n): ", np.log(x),np.log(n)

    m_lum, m_err_lum, b_lum, b_err_lum = linear_fit (np.log(x),np.log(n))
    print "m_lum, m_err_lum: ", m_lum, m_err_lum
    print "b_lum, b_err_lum: ", b_lum, b_err_lum
    print " "

    ax = fig.add_subplot(212)
    #n, bins, patches = ax.hist(f_aper_vec, 45, normed=False, facecolor='yellow', alpha=0.75, label='')
    #plt.xlabel('Fractional luminosity')
    #plt.ylabel('Frequency')
    n, bins, patches = ax.hist(f_auto_vec, 15, normed=True, facecolor='yellow', alpha=0.75, histtype='step', fill=False)
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel('F_AUTO')
    plt.ylabel('Frequency')
    #plt.legend(loc='upper right', fancybox=True, ncol=1, numpoints=1, prop = prop)

    plt.yscale('log',nonposy='clip')
    plt.xscale('log')

    plt.tight_layout()
    pp.savefig()

    np.savetxt ('flux.dat', f_auto_vec)
    np.savetxt ('size.dat', size_vec_kpc)


    size_fit=np.linspace(200, np.max(size_vec_kpc), 1000)
    s_fit= (230000*b_size)*size_fit**m_size

    lum_fit=np.linspace(2000, np.max(f_auto_vec), 1000)
    l_fit= (82*b_lum)*lum_fit**m_lum

    
    fig=plt.figure()
    ax = fig.add_subplot(211)
    ax.errorbar (size_fit, s_fit,yerr=None, fmt='r-', alpha=0.8, label='power-law exponent: \n %g $\pm$ %g' %(m_size, np.sqrt(m_err_size[0])), markersize=3)
    plt.legend(loc='upper left', fontsize=8.5 )
    plt.axvline(x=200, ls='--', c='r')
    n, bins, patches = ax.hist(size_vec_kpc, 22, normed=True, facecolor='green', alpha=0.75, histtype='step', fill=False)
    plt.xlabel('Half-light radius (pc)', size=16)
    #plt.ylabel('Frequency', size=16)
    plt.yscale('log',nonposy='clip')
    plt.xscale('log')
    plt.ylim([1e-4,1e-1])

    

    ax = fig.add_subplot(212)
    ax.errorbar (lum_fit, l_fit,yerr=None, fmt='r-', alpha=0.8, label='power-law exponent: \n %g $\pm$ %g' %(m_lum, np.sqrt(m_err_lum[0])), markersize=3)
    plt.legend(loc='upper right', fontsize =8.5)
    plt.axvline(x=2000, ls='--', c='r')
    n, bins, patches = ax.hist(f_auto_vec, 20, normed=True, facecolor='blue', alpha=0.75, histtype='step', fill=False)
    plt.yscale('log',nonposy='clip')
    plt.xscale('log')
    plt.xlabel('F_AUTO')
    plt.ylim([1e-6,1e-3])
    #plt.ylabel('Frequency')



    #print "size_fit, s_fit: ", size_fit, s_fit
    #print "lum_fit, l_fit: ", lum_fit, l_fit
    plt.tight_layout()


    pp.savefig()


    

    pp.close()




if __name__ == "__main__":
    import pdb, traceback
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
