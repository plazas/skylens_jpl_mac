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


from matplotlib.backends.backend_pdf  import PdfPages


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


output_directory_path="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/OUTPUT_LENSED/"
cmd="ls -d %s*/" %output_directory_path
list_output_dirs=S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
frei_images_path="/Users/amalagon/HST_SL_SIMS/data/frei_galaxy_catalog/"
frei_distances_path="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/z_distances_frei.dat"  # n2403 0.000437 (z) 1931.25 (dA, kpc) 1843.66 (hubble d, kpc)

arcsec_to_rad = 1./206265


def main (argv):
    pp=PdfPages("out.pdf")
    size_vec_pix=[]
    size_vec_kpc=[]
    dA_dict={}
    data=np.genfromtxt(frei_distances_path, dtype=[('gal','|S5'),('z','f8'), ('da','f8'), ('r_hubble','f8')])
    for i, gal in enumerate(data['gal']):
        dA_dict[gal]=np.abs(data['da'][i])
    for dir in list_output_dirs:
        print "dir: ", dir
        print type(dir)
        root=dir.split('/')[-2]
        print "root: ", root
        hd=pf.open(frei_images_path+root+".fits")[0].header
        tel=hd['TELESCOP']
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
        chi2=chi2[4:]  # some weird charcaters at the beginning of the string
        cmd="< %s awk '/host.x/ {print $2}'" %file
        xc=get_number(cmd)
        cmd="< %s awk '/host.y/ {print $2}'" %file
        yc=get_number(cmd)
        cmd="< %s awk '/host.r/ {print $2}'" %file
        r=get_number(cmd)
        cmd="< %s awk '/host.mag/ {print $2}'" %file
        mag=get_number(cmd)
        cmd="< %s awk '/host.n/ {print $2}'" %file
        n=get_number(cmd)
        cmd="< %s awk '/host.q/ {print $2}'" %file
        q=get_number(cmd)
        cmd="< %s awk '/host.pa/ {print $2}'" %file
        beta=get_number(cmd)
        print "model parameters: ", chi2, xc, yc, r, mag, n, q, beta
        # Get model image
        tel=root.split('_')[1][0]
        if tel == 'p':
            pixel_scale=1.19
        else:
            pixel_scale=1.34
        x_size=hd['NAXIS1']
        y_size=hd['NAXIS2']
        x_offset = xc - x_size*1./2
        y_offset = yc - y_size*1./2
        print "offsets: ", x_offset, y_offset
        model_im = draw_sersic_model (mag=mag, r=r, pixel_scale=pixel_scale, q=q, beta=beta, n=n, xc=x_offset, yc=y_offset, x_size=x_size, y_size=y_size)
        # Real data image
        print "real image: ", frei_images_path+root+".fits"
        data=pf.open(frei_images_path+root+".fits")[0].data
        # Difference
        diff=data - model_im
        #Write the diff file to a fits file, and then run Source Extractor
        cmd="rm diff.fits"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        pf.writeto("diff.fits", diff, output_verify='ignore')
        
        output_catalog="OUTPUT_SEXTRACTOR/%s_out.cat"%root
        cmd="mkdir -v OUTPUT_SEXTRACTOR"
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

        cmd="sex diff.fits -c daofind_sex.config -CATALOG_NAME %s" %output_catalog
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        cmd="cp check.fits OUTPUT_SEXTRACTOR/%s_check.fits"%root
        S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
        ## Open sextractor catalog and get the star forming regions
        cat=np.genfromtxt(output_catalog)
        try:
            r=cat[:,8]
        except IndexError:
            continue
        r.sort()
        r=r[r>0.5]
        if len (r) == 0:
            continue
        p=np.percentile (r, 40)  # choose the smaller objects found by Sextractor
        r=r[r<=p]  #radius in pixels
        size_vec_pix.extend(r)
        dA = dA_dict[root[0:5]]
        r_kpc = dA*(r*pixel_scale*arcsec_to_rad)  #dA in kpc
        r_pc=r_kpc*1e3
        size_vec_kpc.extend(r_pc)
        print "Size in pixels; size in pc: ", r, r_pc


        #load the segmentation file
        ap=pf.open("check.fits")[0].data
        
        fig=plt.figure()
        ax=fig.add_subplot (221)
        data, d=histeq(data)
        ax.imshow((data), cmap=cm.gray, origin="lower", interpolation='nearest')
        ax.set_title ("data", size=9)
        ax=fig.add_subplot (222)
        model_im, d=histeq(model_im)
        ax.imshow((model_im), cmap=cm.seismic, origin="lower", interpolation='nearest')
        ax.set_title ("model (chi2: %s)" %chi2, size=9)
        ax=fig.add_subplot (223)
        diff,d=histeq(diff)
        ax.imshow((diff), cmap=cm.gray, origin="lower", interpolation='nearest')
        ax.set_title ("diff", size=9)
        ax=fig.add_subplot (224)
        #ap,d=histeq(ap)
        ax.imshow(ap, cmap=cm.cubehelix , origin="lower", interpolation='nearest', vmin=ap.min()*0.25, vmax=ap.max()*0.3)
        ax.set_title ("apertures of diff", size=9)
        plt.suptitle("%s"%(root))
        plt.tight_layout()

        pp.savefig()

    size_vec_pix=np.array(size_vec_pix)
    size_vec_kpc=np.array(size_vec_kpc)
    print size_vec_pix
    print size_vec_kpc
    size_file=open("star_forming_sizes.dat", 'w')
    for p, x in zip(size_vec_pix, size_vec_kpc):
        line="%g %g \n"%(p,x)
        size_file.write(line)
    size_file.close()

    fig=plt.figure()
    ax = fig.add_subplot(211)
    n, bins, patches = ax.hist(size_vec_pix, 50, normed=False,
                               facecolor='green', alpha=0.75,
                               label='')
    plt.xlabel('Size (pixels)')
    plt.ylabel('histogram')

    ax = fig.add_subplot(212)
    n, bins, patches = ax.hist(size_vec_kpc, 50, normed=False,
                                   facecolor='red', alpha=0.75,
                                   label='')
    plt.xlabel('Size (pc)')
    plt.ylabel('histogram')
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
