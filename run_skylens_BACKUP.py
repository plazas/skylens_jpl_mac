#!/usr/bin/python
import numpy as np
import sys
import galsim

import matplotlib
matplotlib.use('Pdf')

import matplotlib.cm as cm  # color bar, to plot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf  import PdfPages
import matplotlib.patches as patches
import galsim.wfirst as wfirst
from collections import OrderedDict
import astropy.io.fits as pf
import aplpy
import subprocess as S
import os

def run_shell_cmd (cmd):
    print >>sys.stderr, cmd
    S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

path_instruments="/Users/amalagon/share/astro/Instruments/"
path_filters="/Users/amalagon/share/astro/Filters/"
path_psf="/Users/amalagon/share/astro/PSF/"
path_skylens="/Users/amalagon/git/skylens3/"
path_source_cat=path_skylens+"python/"

class instrument:
    def __init__(self, name, pixel_scale, dict_all, fov, rgb):
        self.name = name
        self.pixel_scale = pixel_scale
        self.dict_all=dict_all
        self.fov=fov
        self.rgb=rgb

###############
name='JWST_NIRCAM_SHORT_LAMBDA'
dict_all={'F070W':{'tp':path_filters+"F070W_NRC_and_OTE_ModAB_mean_angstroms_JWST.ecf",'psf':path_psf+"JWST_PSF_F070W.fits", 'cat':path_source_cat+"mycat_JWST_NIRCam_F070W.cat", 'sky': 0.20, 'time': 500000 , 'seeing': 0.0}, \
          'F090W':{'tp': path_filters+"F090W_NRC_and_OTE_ModAB_mean_angstroms_JWST.ecf",'psf': path_psf+"JWST_PSF_F090W.fits",'cat':path_source_cat+"mycat_JWST_NIRCam_F090W.cat",'sky': 0.318, 'time': 500000 , 'seeing': 0.0} , \
          'F115W':{'tp': path_filters+"F115W_NRC_and_OTE_ModAB_mean_angstroms_JWST.ecf",'psf': path_psf+"JWST_PSF_F115W.fits",'cat':path_source_cat+"mycat_JWST_NIRCam_F115W.cat",'sky': 0.315, 'time': 500000 , 'seeing': 0.0} , \
          'F150W':{'tp': path_filters+"F150W_NRC_and_OTE_ModAB_mean_angstroms_JWST.ecf",'psf': path_psf+"JWST_PSF_F150W.fits",'cat':path_source_cat+"mycat_JWST_NIRCam_F150W.cat",'sky': 0.321, 'time': 500000 , 'seeing': 0.0} , \
          'F200W':{'tp': path_filters+"F200W_NRC_and_OTE_ModAB_mean_angstroms_JWST.ecf",'psf': path_psf+"JWST_PSF_F200W.fits",'cat':path_source_cat+"mycat_JWST_NIRCam_F200W.cat",'sky': 0.289, 'time': 500000 , 'seeing': 0.0} }
ps=0.032/4  #scale of the oversampled HDU of the PSF file.
fov= 204 #arcsec
rgb=['F200W','F150W','F090W']
jw = instrument (name, ps, dict_all,fov,rgb)

################
#### H158 is red, J129 is green, and Y106 is blue. 
name='WFIRST_WFI'
dict_all={'Y106':{'tp':path_filters+"WFIRST_Y106.ecf",'psf': path_psf+"WFIRST_Y106_PSF_from_galsim.fits",'cat': path_source_cat+"mycat_WFIRST_Y106.cat",'sky': 0.6547, 'time': 500000 , 'seeing': 0.0}, \
          'J129':{'tp':path_filters+"WFIRST_J129.ecf",'psf': path_psf+"WFIRST_J129_PSF_from_galsim.fits",'cat':path_source_cat+"mycat_WFIRST_J129.cat",'sky': 0.6696,'time': 500000 , 'seeing': 0.0} , \
          'H158':{'tp':path_filters+"WFIRST_H158.ecf",'psf': path_psf+"WFIRST_H158_PSF_from_galsim.fits",'cat':path_source_cat+"mycat_WFIRST_H158.cat",'sky': 0.6543,'time': 500000 , 'seeing': 0.0} , \
          'F184':{'tp':path_filters+"WFIRST_F184.ecf",'psf': path_psf+"WFIRST_F184_PSF_from_galsim.fits",'cat':path_source_cat+"mycat_WFIRST_F184.cat",'sky': 1.51292,'time': 500000 , 'seeing': 0.0} }
ps=0.11/3
fov= 204 #arcsec
rgb=['Y106','J129','H158']
wf = instrument (name, ps, dict_all, fov, rgb)

################
name='HST_ACS_WFC'
dict_all={'F435W':{'tp':path_filters+"F435W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F435W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F435W.cat",'sky': 0.0364,'time': 500000 , 'seeing': 0.0}, \
          'F475W':{'tp':path_filters+"F475W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F475W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F475W.cat", 'sky':0.0621,'time': 500000 , 'seeing': 0.0} , \
          'F606W':{'tp':path_filters+"F606W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F606W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F606W.cat",'sky': 0.1325,'time': 500000 , 'seeing': 0.0} , \
          'F625W':{'tp':path_filters+"F625W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F625W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F625W.cat",'sky': 0.0877,'time': 500000 , 'seeing': 0.0}, \
          'F775W':{'tp':path_filters+"F775W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F775W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F775W.cat",'sky': 0.0831,'time': 500000 , 'seeing': 0.0} , \
          'F814W':{'tp':path_filters+"F814W_WFC.res",'psf': path_psf+"HST_ACS_WFC_F814W_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F814W.cat",'sky': 0.108,'time': 500000 , 'seeing': 0.0} , \
          'F850LP':{'tp':path_filters+"F850LP_WFC.res",'psf': path_psf+"HST_ACS_WFC_F850LP_PSF00.fits",'cat':path_source_cat+"mycat_HST_WFC_F850LP.cat",'sky': 0.0453,'time': 500000 , 'seeing': 0.0}}
ps= 0.0495
fov= 204 #arcsec
rgb=['F775W','F625W','F475W']
hst = instrument (name, ps, dict_all, fov,rgb)


############## LUVOIR
#### Using Optical PSF, and using HST throughput curves instead. Setting sky to zero, so skylens calculates it.
#### Do 3 wavelengths: 475 (blue), 625 (green), 775 (red) for the optical PSF and the WFC filters. Assume 0.0 obscuration, 15.1m telescope diameter, no aberrations...Airy, basically.
name='LUVOIR_HDI'
dict_all={
          'F475W':{'tp':path_filters+"F475W_WFC.res",'psf': path_psf+"LUVOIR_HDI_PSF_GALSIM_OPT_475.fits",'cat':path_source_cat+"mycat_HST_WFC_F475W.cat", 'sky':0.0,'time': 500000 , 'seeing': 0.0} , \
          
          'F625W':{'tp':path_filters+"F625W_WFC.res",'psf': path_psf+"LUVOIR_HDI_PSF_GALSIM_OPT_625.fits",'cat':path_source_cat+"mycat_HST_WFC_F625W.cat",'sky': 0.0,'time': 500000 , 'seeing': 0.0}, \
          'F775W':{'tp':path_filters+"F775W_WFC.res",'psf': path_psf+"LUVOIR_HDI_PSF_GALSIM_OPT_775.fits",'cat':path_source_cat+"mycat_HST_WFC_F775W.cat",'sky': 0.0,'time': 500000 , 'seeing': 0.0} }
ps= 0.00274 
fov= 204 #arcsec
rgb=['F775W','F625W','F475W']
lv = instrument (name, ps, dict_all, fov,rgb)


################
## PSF from /Users/amalagon/share/astro/PSF/HSC_PSF_WIDE_AIHARA17_RA180DEG_DEC0DEG. I renamed them to HSC_SUBARU_g.fits etc
## Seeing from Aihara et al 2017: 0.72 0.67 0.56 0.63 0.64
name='SUBARU_HSC'
dict_all={'g':{'tp':path_filters+"hsc_g.ecf",'psf': path_psf+"HSC_SUBARU_PSF_g.fits",'cat':path_source_cat+"mycat_HSC_SUBARU_g.cat", 'sky':35.08, 'time': 500000 , 'seeing': 0.72}, \
          'r':{'tp':path_filters+"hsc_r.ecf",'psf': path_psf+"HSC_SUBARU_PSF_r.fits",'cat':path_source_cat+"mycat_HSC_SUBARU_r.cat", 'sky':48.54, 'time': 500000 , 'seeing': 0.67}, \
          'i':{'tp':path_filters+"hsc_i.ecf",'psf': path_psf+"HSC_SUBARU_PSF_i.fits",'cat':path_source_cat+"mycat_HSC_SUBARU_i.cat", 'sky':75.74, 'time': 500000 , 'seeing': 0.56}, \
          'z':{'tp':path_filters+"hsc_z.ecf",'psf': path_psf+"HSC_SUBARU_PSF_z.fits",'cat':path_source_cat+"mycat_HSC_SUBARU_z.cat", 'sky':45.60, 'time': 500000 , 'seeing': 0.63}, \
          'y':{'tp':path_filters+"hsc_y.ecf",'psf': path_psf+"HSC_SUBARU_PSF_y.fits",'cat':path_source_cat+"mycat_HSC_SUBARU_y.cat", 'sky':100.15,'time': 500000 , 'seeing': 0.64} }

ps=0.17
fov=1800
rgb=['z','i','r']
hsc = instrument (name, ps, dict_all, fov,rgb)


###############
### use background from exposure time calculator for DECam
### For seeing, use median PSF FWHM from DR1 of DES : 1.119 0.958 0.880 0.836 0.904
name='LSST'
dict_all={'g':{'tp':path_filters+"LSST_g.res", 'psf':  path_psf+"LSST_PSF_GALSIM_OPT_KOLMO_g_500.fits",'cat':path_source_cat+"mycat_LSST_g.cat",'sky': 60.26, 'time': 500000 , 'seeing':  1.119}, \
          'r':{'tp':path_filters+"LSST_r.res", 'psf':  path_psf+"LSST_PSF_GALSIM_OPT_KOLMO_r_650.fits",'cat':path_source_cat+"mycat_LSST_r.cat",'sky': 22.79, 'time': 500000 , 'seeing':  0.958}, \
          'i':{'tp':path_filters+"LSST_i.res", 'psf':  path_psf+"LSST_PSF_GALSIM_OPT_KOLMO_i_750.fits",'cat':path_source_cat+"mycat_LSST_i.cat",'sky': 41.57, 'time': 500000 , 'seeing':  0.880}, \
          'z':{'tp':path_filters+"LSST_z.res", 'psf': path_psf+"LSST_PSF_GALSIM_OPT_KOLMO_z_860.fits",'cat':path_source_cat+"mycat_LSST_r.cat",'sky': 116.75,'time': 500000 , 'seeing':  0.836}, \
          'y':{'tp':path_filters+"LSST_y.res", 'psf': path_psf+"LSST_PSF_GALSIM_OPT_KOLMO_y_1000.fits",'cat':path_source_cat+"mycat_LSST_y.cat",'sky': 35.12, 'time': 500000 , 'seeing': 0.904} }
ps=0.2
fov=1800
rgb=['z','i','r']
lsst = instrument (name, ps, dict_all, fov,rgb)



######### Now run sky lens for the configurations above 

to_run=[jw, wf] #, hst, lv, hsc, lsst]
directory_suffix="TEST1"

for i in to_run:
    ### create directory where to run skylens and store images
    dire="%s_%s"%(i.name, directory_suffix)
    cmd="mkdir -v %s"%(dire); run_shell_cmd(cmd)
    #cmd="cd %s"%dire; run_shell_cmd(cmd)
    
    for band in i.dict_all:
        print "RUNNING: ", i.name, band


        ###modify configuration file of skylens  Use: cat file | sed 's/^.*STRING_TO_BE_REPLACED.*$/LINE_TO_REPLACE_WITH/' >new_file
        new_skylens_file="skylens_andres_%s_%s.inp" %(i.name, band)
        cmd="cp skylens_andres_template.inp %s/%s" %(dire,new_skylens_file); run_shell_cmd(cmd)
        
        #old_string, new_string = "catalog", "%s"i.dict_all[band]['cat']
        #cmd="< %s/%s sed -i 's/^.*%s.*$/%s/'" %(dire, new_skylens_file, old_string, new_string) ; run_shell_cmd(cmd)
        """
        old_string, new_string = "size of the field", "%s" %i.fov
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)   
        
        old_string, new_string = "name of the telescope", "%s" %i.name
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)   

        old_string, new_string = "exposure time", "%s" %i.dict_all[band]['time']
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)

        old_string, new_string = "seeing", "%s" %i.dict_all[band]['seeing']
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)

        old_string, new_string = "psf image", "%s" %i.dict_all[band]['psf']
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)

        old_string, new_string = "psf scale", "%s" %i.pixel_scale
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)

        old_string, new_string = "throughput", "%s" %i.dict_all[band]['tp']
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %(old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)

        old_string, new_string = "sky surf brightness", "%s" %i.dict_all[band]['sky']
        cmd="sed -i '' 's/^.*%s.*$/%s/' %s/%s" %( old_string, new_string, dire,new_skylens_file) ; run_shell_cmd(cmd)
        """
        #### RUN SKYLENS HERE
        cmd="time "+path_skylens+"SkyLens_v3.0_bkg_andres.x -i %s/%s"%(dire,new_skylens_file) #; run_shell_cmd(cmd)
        #print "CMD: ", cmd
        #args=['/Users/amalagon/git/skylens3/SkyLens_v3.0_bkg_andres.x', '-i=skylens_andres_JWST_NIRCAM_SW_F200W.inp']
        #S.call(args) 
        os.system("/Users/amalagon/git/skylens3/SkyLens_v3.0_bkg_andres.x -i JWST_NIRCAM_SHORT_LAMBDA_TEST1/skylens_andres_JWST_NIRCAM_SHORT_LAMBDA_F200W.inp")
        #S.Popen(["/Users/amalagon/git/skylens3/SkyLens_v3.0_bkg_andres.x"], shell=True, stdout=S.PIPE).communicate("-i skylens_andres_JWST_NIRCAM_SW_F200W.inp")
        sys.exit()
        #### Rename output from skylens: mult_images.cat, mult_images.reg (if aplicable), testimage.fits.dat, testimage.fits, testimage_noisy.fits.dat, testimage_noisy.fits, journal.log
        cmd="mv %s/mult_images.cat %s/mult_images_%s_%s.cat"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/mult_images.reg %s/mult_images_%s_%s.reg"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/testimage.fits.dat %s/testimage_%s_%s.fits.dat"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/testimage.fits %s/testimage_%s_%s.fits"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/testimage_noisy.fits %s/testimage_noisy_%s_%s.fits"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/testimage_noisy.fits.dat %s/testimage_noisy_%s_%s.fits.dat"%(dire,dire,i.name, band); run_shell_cmd(cmd)
        cmd="mv %s/journal.log %s/journal_%s_%s.log"%(dire,dire,i.name, band); run_shell_cmd(cmd)

    ####Go back to parent directory
    #cmd="cd .."; run_shell_cmd(cmd)
        

#### Now, use 3 images ("RGB") from each run and combine them in a single color image. 
def combine_rgb (instrument,dire):
    #noisy
    file_r=dire+"testimage_noisy_%s_%s.fits"%(instrument.name,instrument.rgb[0])
    file_g=dire+"testimage_noisy_%s_%s.fits"%(instrument.name,instrument.rgb[1])
    file_b=dire+"testimage_noisy_%s_%s.fits"%(instrument.name,instrument.rgb[2])
    cutout_r=pyfits.open(file_r)[0].data
    cutout_g=pyfits.open(file_g)[0].data
    cutout_b=pyfits.open(file_b)[0].data

    image_cube = np.zeros((3, cutout_r.shape[0], cutout_r.shape[1]), dtype=np.float32)
    image_cube[0,:,:]=cutout_r
    image_cube[1,:,:]=cutout_g
    image_cube[2,:,:]=cutout_b

    pyfits.writeto(dire+'new.fits', image_cube, clobber=True)

    aplpy.make_rgb_image(dire+'new.fits',dire+"combined_rgb_noisy_%s_%s_%s_%s.png"%(instrument.name,instrument.rgb[0],instrument.rgb[1],instrument.rgb[2]))

    # no noise 
    file_r=dire+"testimage_%s_%s.fits"%(instrument.name,instrument.rgb[0])
    file_g=dire+"testimage_%s_%s.fits"%(instrument.name,instrument.rgb[1])
    file_b=dire+"testimage_%s_%s.fits"%(instrument.name,instrument.rgb[2])
    cutout_r=pyfits.open(file_r)[0].data
    cutout_g=pyfits.open(file_g)[0].data
    cutout_b=pyfits.open(file_b)[0].data

    image_cube = np.zeros((3, cutout_r.shape[0], cutout_r.shape[1]), dtype=np.float32)
    image_cube[0,:,:]=cutout_r
    image_cube[1,:,:]=cutout_g
    image_cube[2,:,:]=cutout_b

    pyfits.writeto(dire+'new.fits', image_cube, clobber=True)

    aplpy.make_rgb_image(dire+'new.fits',dire+"combined_rgb_%s_%s_%s_%s.png"%(instrument.name,instrument.rgb[0],instrument.rgb[1],instrument.rgb[2]))
    

print "Making a composite image of 3 colors 'RGB'."
for i in to_run:
    dire="%s_%s"%(i.name, directory_suffix)
    combine_rgb(i, dire)


print "End of 'run_skylens.py'"    
