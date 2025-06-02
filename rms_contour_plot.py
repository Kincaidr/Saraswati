import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS

def plot(fits_image, path, name):

    hdu = fits.open(fits_image)[0]
    header = fits.getheader(fits_image)
    wcs = WCS(header)
    wcs=wcs[0,0,:,:]

    image_data = hdu.data

    if len(image_data.shape) > 2:
        image_data = image_data[0, 0,:, :] 
        
    smoothed_data = gaussian_filter(image_data, sigma=4)*1e6
    sigma=11
    levels=[sigma*i for i in [2,4,8,16,32,64]]
    vmin=1*sigma
    vmax=10*sigma

    fig, ax = plt.subplots(figsize=(8, 8),subplot_kw={'projection': wcs})
    cax = ax.imshow(smoothed_data, cmap='viridis_r', origin='lower', interpolation='none', vmin=vmin, vmax=vmax)
    contour = ax.contour(smoothed_data, levels=levels, colors='black')

    colorbar=fig.colorbar(cax, ax=ax, orientation='vertical', label=r'RMS ($\mu$Jy/Beam)')
    colorbar.set_label(r'RMS ($\mu$Jy/Beam)', fontsize=18, color='black')
    ax.set_xlabel('RA (J2000)',size=17 )
    ax.set_ylabel('Dec (J2000)', size=17)
    ax.set_xlim(1200, 4800)  # Set x-axis limits (change to your desired range)
    ax.set_ylim(1200, 4800)  #
    plt.savefig(outname+name+'_rms_contour_map.png', bbox_inches='tight', pad_inches=0)
    plt.show() 

if "__main__":
    name='Zwcl2341'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    outname='/home/kincaid/Desktop/Saraswati_codes/plots/'
    fits_image = path+name+'_full_rms_map.fits'

    plot(fits_image, path, name)