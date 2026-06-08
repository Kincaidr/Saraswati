
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from astropy.cosmology import Planck18 as cosmo
from astropy.wcs.utils import pixel_to_skycoord
from matplotlib.patches import Circle
from regions import Regions



def scale_bar(z,pix_scale):
    D_A = cosmo.angular_diameter_distance(z).value
    scale_2Mpc_arcsec = (2 / D_A) * (180 / 3.141592) * 3600 
    scale= scale_2Mpc_arcsec / pix_scale
    return(scale)

def image_analyse(fits_image):
    hdu = fits.open(fits_image)[0]
    header = fits.getheader(fits_image)
    wcs = WCS(header)
    wcs=wcs[0,0,:,:]
 
    image_data = hdu.data
    if len(image_data.shape) > 2:
        image_data = image_data[0, 0,:, :] *1e6
        for key in list(header.keys()):
            if '3' in key or '4' in key:
                header.remove(key, ignore_missing=True, remove_all=True)

        header['NAXIS'] = 2
        header['NAXIS1'] = image_data.shape[1]
        header['NAXIS2'] = image_data.shape[0]

    # sigma=10
    # vmin=1*sigma
    # vmax=20*sigma
    fits_image_new= fits_image.replace('.fits', '_new.fits')
    fits.writeto(fits_image_new, image_data, header, overwrite=True)
    return(fits_image_new,header)

# def plot(fits_image, plots, outname):

#     image_data,vmin, vmax, wcs=image_analyse(fits_image)
#     fig, ax = plt.subplots(figsize=(8, 8))
#     #ax = fig.add_subplot(1, 1, 1, projection=wcs)  
#     cax1 = ax.imshow(image_data * 1e6, cmap='gray', origin='lower', interpolation='none', vmin=vmin, vmax=vmax)
#     ax.set_xlabel('RA (J2000)', size=16)
#     ax.set_ylabel('Dec (J2000)', size=16)
#     colorbar=fig.colorbar(cax1, orientation='vertical')
#     colorbar.set_label(r'($\mu$Jy/Beam)', fontsize=16, color='black')
#     ax.set_xlim(1200, 4800)  # Set x-axis limits (change to your desired range)
#     ax.set_ylim(1200, 4800)  #
#     plt.tight_layout(rect=[0, 0, 0.95, 1]) 
#     plt.tick_params(axis='both', which='major', labelsize=13, length5, width=1)  # Increase size of major tick labels
#     plt.tick_params(axis='both', which='minor', labelsize=13, length=5, width=1) # Adjust layout to fit colorbar
#     plt.savefig(plots+outname+'_full_image.png')
#     plt.show()
    
import aplpy
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

def plot(fits_image1):
    # Load and analyse the image
    new_fits,header= image_analyse(fits_image1)
    ra_center, dec_center = header['CRVAL1'],header['CRVAL2'] # Get the center coordinates from the WCS
    fig = plt.figure(figsize=(8, 8))
    f = aplpy.FITSFigure(new_fits, figure=fig)
    f.show_colorscale(cmap='inferno_r', pmin=1, pmax=99.5,)
    f.recenter(ra_center, dec_center,radius=0.4)
    f.axis_labels.set_font(size=16)
    f.tick_labels.set_font(size=13)
    f.show_regions(regfile)
    f.axis_labels.set_xtext('RA (J2000)')
    f.axis_labels.set_ytext('Dec (J2000)')
    # Add a colorbar
    f.add_colorbar()
    f.colorbar.set_axis_label_text(r'($\mu$Jy/Beam)')
    f.colorbar.set_axis_label_font(size=16)
    f.colorbar.set_font(size=13)

    # Add a scale bar
    # f.add_scalebar((scale_bar(z, pix_scale)))
    # f.scalebar.set_label('500 kpc')  # Adjust as needed
    # f.scalebar.set_color('white')
    # f.scalebar.set_linewidth(2)
    # f.scalebar.set_font_size(12)

    # Save the plot
    output_file = '/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/VLA_image.png'
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()


if "__main__":
    z=0.27
    pix_scale=1.5
    regfile='A2631_VLA_circle_regions.reg'
    name='A2631'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    fits_image1 =  'oussid._Abell_2631__sci.L_band.cont.selfcal.I.pbcor.tt0.fits.img.conv.fits'
    plot(fits_image1)
