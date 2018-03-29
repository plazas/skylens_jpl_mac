#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import esutil
from scipy import ndimage
import sys
import matplotlib.cm as cm
import galsim
import subprocess as S
import astropy.io.fits as pf
import os



## Input: directory with FITS files of Frei et al. galaxies.
## Description: Open each file, run LENSED to create Sersic models, using values in header as input for priors.
## Output: standard output of LENSED (in a directory), plus LENSED requested output.

## - 3 <= T <= 9

# Header example
"""
    SIMPLE  =                    T
    BITPIX  =                   16
    NAXIS   =                    2
    NAXIS1  =                  721
    NAXIS2  =                  721
    OBJECT  = 'NGC 2403  '
    FILTER  = 'GC1                    (4805/75)'/ 3 inch filter
    RA      = '07:32:02.60'
    DEC     = '+65:42:44.2'
    AIRMASS =         1.462037E+00
    HISTORY   'Flat-fielded, star-subtracted, calibrated, header-corrected '
    TELESCOP= 'Palomar 1.5'
    INSTRUME= 'Wide field'
    OBSERVER= 'J. E. Gunn'
    AUTHOR  = 'Z. Frei et al. '
    REFERENC= 'AJ, 1995, in press '
    DATE    = '12/07/95  '
    EQUINOX =         1.950000E+03
    BUNIT   = 'COUNTS    '
    HA      = '-03:52:07.9'                / Hour Angle
    LST     = ' 11:28:10.3'                / Sidereal Time
    TIME    = '04:25:24.2'
    FILENUM =                   16
    DECORDER=                    F
    BSCALE  =         1.000000E+00
    BZERO   =         0.000000E+00
    CRVAL1  =         1.130083E+02
    CRVAL2  =         6.571222E+01
    CRPIX1  =         3.680000E+02
    CRPIX2  =         3.560000E+02
    CDELT1  =        -3.305555E-04
    CDELT2  =         3.305555E-04
    CTYPE1  = 'RA---TAN  '
    CTYPE2  = 'DEC--TAN  '
    DATAMAX =         1.281400E+04
    DATAMIN =         1.400000E+02
    CCDNUM  =                    1
    EXPOSURE=         6.000000E+01
    FILENAME= 'n2403_pr.fits'
    FILTER1 = 'r         '
    SKY     =         1.556470E+03
    SKYSIG  =         4.045702E+01
    XOFFSET =                   54
    YOFFSET =                   58
    G_CENT_X=                  367
    G_CENT_Y=                  355
    SATURAT =                    1
    DNAT0_ST=         1.382728E+11
    DNAT0_BV=         1.334337E+11
    B_RC3   =         8.930000E+00
    B-V_RC3 =         4.700000E-01
    PSF_FWHM=         2.288388E+00
    PA_RC3  =         1.270000E+02
    BOA_RC3 =         2.500000E-01
    SIZE_RC3=         2.340000E+00
    VELO_RC3=         1.299699E+02                                                  
    TYPE_RC3= '.SXS6..   '                                                          
    DATE-OBS= '05/05/91  '
"""

import ConfigParser



directory_path="/Users/amalagon/HST_SL_SIMS/data/frei_galaxy_catalog/"
template_ini_path="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/"
def run_shell_cmd (cmd):
    print >>sys.stderr, cmd
    #S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()
    S.Popen([cmd], shell=True, stdout=S.PIPE).stdout.read()

general_output_dir_path="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/OUTPUT_LENSED_BULGE_PLUS_DISK/"
cmd="mkdir -v %s" %general_output_dir_path
run_shell_cmd (cmd)


if not os.path.isdir(general_output_dir_path):
    os.mkdir(general_output_dir_path)

#file_name = os.path.join('output','demo1_gal.fits')


cmd="ls %s*.fits" %directory_path
list=S.Popen([cmd], shell=True, stdout=S.PIPE).communicate()[0].split()

LENSED_config_files_dir_path="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/CONFIG_FILES_LENSED_BULGE_PLUS_DISK/"
data_path="/Users/amalagon/HST_SL_SIMS/data/frei_galaxy_catalog/"
LENSED_executable="/Users/amalagon/git/lensed/bin/lensed"
model=2  #Either 1 (single Sersic) or 2 (bulge+disk, both Sersic)


def main (argv):
    for file in list:
        print "File: ", file
        hd=pf.open(file)[0].header
        tel=hd['TELESCOP']
        if tel == 'Palomar 1.5':
            pixel_scale = 1.19  #arcsec per pixel
        elif tel == 'Lowell 1.1':
            pixel_scale = 1.35
        else:
            print "Warning: Telescope is not Palomar nor Lewis. Setting plate scale to 1."
            pixel_scale=1.0
        
        time=hd['EXPOSURE']
        sky=hd['SKY']
        x0=hd['G_CENT_X']
        y0=hd['G_CENT_Y']
        pa0=hd['PA_RC3']  # From north to east (up to left on ds9)
        if pa0 <= 90:   # I think pa in LENSED is from x-axis
            pa0+=90
        else:
            pa0-=90
        r0=hd['SIZE_RC3']  # log 10 of size estimate
        r0=10**r0*1./pixel_scale  # arcseconds --> pixels

        root=file.split('/')[-1].split('.')[0]
        #Now modify the .ini file for this particular image
        Config = ConfigParser.ConfigParser()
        if model == 1:
            Config.read(template_ini_path+"template.ini")
        elif model == 2:
            Config.read(template_ini_path+"template-bulge-disk.ini")
        else:
            print "ERROR in template: 'model' is either 1 (single Sersic) or '2' (bulge+disk, both Sersic). "
            sys.exit(1)

        Config.set('options','image', data_path+root+".fits")
        Config.set('options','gain', 2)
        Config.set('options','output', "true")
        particular_output_dir=root+"/"
        if not os.path.isdir(general_output_dir_path+particular_output_dir):
            os.mkdir(general_output_dir_path+particular_output_dir)
        Config.set('options','root', general_output_dir_path+particular_output_dir+root)
        Config.set('objects','sky', "sky")
        if model == 1:
            Config.set('objects','host', "sersic")
            Config.set('priors','host.x', "image unif %g %g" %(x0-5, x0 + 5))
            Config.set('priors','host.y', "image unif %g %g" %(y0-5, y0 + 5))
            Config.set('priors','host.r', "unif %g %g" %(r0-10 , r0 + 10))
            Config.set('priors','host.mag', "unif %g %g" %(-22, -10))
            Config.set('priors','host.n', "unif %g %g" %(0.5, 4))
            Config.set('priors','host.q', "unif %g %g" %(0.3, 0.9))
            Config.set('priors','host.pa', "wrap unif %g %g" %(pa0 - 10, pa0 + 10))
        elif model == 2:
            Config.set('objects','bulge', "sersic")
            Config.set('objects','disk', "sersic")

            Config.set('priors','bulge.x', "image unif %g %g" %(x0-5, x0 + 5))
            Config.set('priors','bulge.y', "image unif %g %g" %(y0-5, y0 + 5))
            Config.set('priors','bulge.r', "unif %g %g" %(r0-10 , r0 + 10))
            Config.set('priors','bulge.mag', "unif %g %g" %(-22, -10))
            Config.set('priors','bulge.n', "unif %g %g" %(0.5, 4))
            Config.set('priors','bulge.q', "unif %g %g" %(0.3, 0.9))
            Config.set('priors','bulge.pa', "wrap unif %g %g" %(pa0 - 10, pa0 + 10))

            Config.set('priors','disk.x', "image unif %g %g" %(x0-5, x0 + 5))
            Config.set('priors','disk.y', "image unif %g %g" %(y0-5, y0 + 5))
            Config.set('priors','disk.r', "unif %g %g" %(r0-10 , r0 + 10))
            Config.set('priors','disk.mag', "unif %g %g" %(-22, -10))
            Config.set('priors','disk.n', "unif %g %g" %(0.5, 4))
            Config.set('priors','disk.q', "unif %g %g" %(0.3, 0.9))
            Config.set('priors','disk.pa', "wrap unif %g %g" %(pa0 - 10, pa0 + 10))
        else:
            print "ERROR: 'model' is either 1 (single Sersic) or '2' (bulge+disk, both Sersic). "
            sys.exit(1)
        
        Config.set('priors','sky.bg', sky)
        
        cfgfile = open(root+".ini",'w')
        Config.write(cfgfile)
        cfgfile.close()
        if not os.path.isdir(LENSED_config_files_dir_path):
            os.mkdir(LENSED_config_files_dir_path)
        cmd="mv %s.ini %s" %(root, LENSED_config_files_dir_path)
        run_shell_cmd(cmd)

        #run lensed
        print " "
        print "Before running LENSED for %s.fits " %root
        cmd="%s %s.ini | tee %s%s-stdot.dat" %(LENSED_executable, LENSED_config_files_dir_path+root,general_output_dir_path+particular_output_dir, root)
        run_shell_cmd(cmd)
        print " "


if __name__ == "__main__":
    import pdb, traceback
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)




