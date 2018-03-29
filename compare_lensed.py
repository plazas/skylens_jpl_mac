#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import esutil
from scipy import ndimage
import sys
import galsim
import matplotlib.cm as cm
import galsim



def get_magnitude_lensed (f_tot):
    return np.log10 (f_tot)/(-0.4)


def get_flux_lensed (mag):
    return 10**(-0.4*mag)


def main(argv):
    # Read in real galaxy
    import fitsio
    file_name="/Users/amalagon/HST_SL_SIMS/FREI_GUNN_GALAXIES_LENSED/n4303_pr.fits"
    data=fitsio.read (file_name)
    print "data.shape", data.shape

    base_size_x, base_size_y = data.shape[0], data.shape[1]

    # Create model with parameters calculated by LENSED

    """
    Dim No.       Mean        Sigma
    1    0.155640000000000214E+04                         NaN
    2    0.000000000000000000E+00    0.000000000000000000E+00
    3    0.000000000000000000E+00    0.000000000000000000E+00
    4    0.361035705890123097E+03    0.173175693747794103E-01
    5    0.355279701709864412E+03    0.146791022921851227E-01
    6    0.999998468534430458E+02    0.152707053384766361E-03
    7   -0.193274297910504060E+02    0.223224607045188018E-03
    8    0.110942613348518315E+01    0.517501148686420887E-03
    9    0.602248424972319785E+00    0.251145232822460945E-03
    10   0.332000808207795259E+02    0.214766882134241575E-01
    """



    flux_gal=get_flux_lensed(-17.6176)    # Lensed documentation: L_tot = 10**(-0.4m)
    print "flux: ", flux_gal
    #flux_gal*=100000000
    size_gal=54.8553   # In arcseconds
    pixel_scale=1.19
    sky_level =594.4662
    q_gal=0.8387
    beta_gal= 106.7037
    n=1.4175
    random_seed=1534225
    rng = galsim.BaseDeviate(random_seed)
    
    gal=galsim.Sersic(n=n, half_light_radius=size_gal, flux=flux_gal).shear(q=q_gal, beta=beta_gal*galsim.degrees)
    model_lensed_im=gal.drawImage(image=galsim.ImageF(base_size_x, base_size_y),scale=pixel_scale)
    # Add Poisson noise to the image:
    model_lensed_im.addNoise( galsim.PoissonNoise(rng, sky_level * pixel_scale**2))


    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(14, 14))
    ax1.imshow(model_lensed_im.array, cmap=cm.seismic)
    ax1.set_title ("Model Galaxy - Galsim")
    ax2.imshow(np.arcsinh(data), cmap=cm.seismic)
    ax2.set_title ("Actual Galaxy - Parameters fitted by Lensed")
    diff = data - model_lensed_im.array
    fitsio.write ("diff.fits", diff)
    ax3.imshow(np.arcsinh (diff), cmap=cm.seismic)
    ax3.set_title ("difference")
    """
    
    #fig.savefig("out.png")
    #fig.show()
    
    diff = data - model_lensed_im.array
    fig, (ax) = plt.subplots (nrows=1, ncols=1, figsize=(7,7))
    ax.imshow(np.arcsinh (diff), cmap=cm.seismic)
    fig.savefig("diff.png")
    stop

    #cax3=ax3.imshow( (exp_lensed_im.array - exp_sim_im.array)/exp_sim_im.array, cmap=cm.seismic)
    #ax3.set_title ("Fractional Difference")
    #fig.colorbar(cax3)






    ##### Now let's try a real galaxy from COSMOS postage stamps that come with galSim

    ## Read in real COSMOS postage stamp
    import fitsio
    file_name="/Users/amalagon/HST_SL_SIMS/lensed-tests/cosmos0.fits"
    data=fitsio.read (file_name)
    print "data.shape", data.shape



    ### Results by lensed (using cosmos29.ini)


    #parameter           mean       sigma          ML         MAP
    #------------------------------------------------------------
    #    host.x           29.9328      2.8799     27.6539     29.0682
    #    host.y           59.8873      2.8952     64.4414     55.8267
    #    host.rs           8.9723      1.1239      8.2299      8.8282
    #    host.mag          6.2480      7.2542     -4.0550     19.0027
    #    host.q            0.7080      0.1154      0.5692      0.7231
    #    host.pa         135.2280     26.1256    109.5517    155.1082



    # Create Exponential model with parameters from above
    
    flux_gal=get_flux_lensed(8.2983)    # Lensed documentation: L_tot = 10**(-0.4m)
    size_gal=8.1964
    pixel_scale=1.0
    sky_level =10.0
    q_gal=0.4707
    beta_gal=-2
    #base_size = 120

    #image=galsim.ImageF(data.shape[0], data.shape[1])

    big_fft_params = galsim.GSParams(maximum_fft_size=4096*6)
    #galsim.GSParams.maximum_fft_size = 10000
    gal= galsim.Convolve(  galsim.Exponential(flux=flux_gal, scale_radius=size_gal).shear(q=q_gal, beta=beta_gal*galsim.degrees) , galsim.Pixel(pixel_scale), gsparams=big_fft_params )
    exp_cosmos_lensed_im=gal.drawImage(image=galsim.ImageF(data.shape[1], data.shape[0]), scale=pixel_scale, method='no_pixel')
    #exp_cosmos_lensed_im.addNoise (galsim.getCOSMOSNoise (file_name=file_name, rng=rng))
    
    


    fig, (ax5, ax6) = plt.subplots (nrows=1, ncols=2, figsize = (14, 7))
    ax6.imshow (np.arcsinh(exp_cosmos_lensed_im.array), cmap=cm.seismic)
    ax6.set_title ("COSMOS gal: Fitted by lensed (exponential), drawn in GalSim")
    ax5.imshow (np.arcsinh(data), cmap=cm.seismic)
    ax5.set_title ("COSMOS gal" )
    fig.savefig("exponential_cosmos29_vs_lensed.png")
    fig.show()
    stop






if __name__ == "__main__":
    import pdb, traceback
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)



