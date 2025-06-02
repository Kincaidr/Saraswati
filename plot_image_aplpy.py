
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from astropy.cosmology import Planck18 as cosmo
import aplpy
import os.path as op

def scale_bar(z,pix_scale):
    D_A = cosmo.angular_diameter_distance(z).value
    scale_2Mpc_arcsec = (2 / D_A) * (180 / 3.141592) * 3600 
    scale= scale_2Mpc_arcsec / pix_scale
    return(scale)   


def get_fitsimage(fits_image,new_fits_image):
    with fits.open(fits_image) as hdul:
        image_data = hdul[0].data
        image_data=image_data[0,0,:,:]*1e6
        header=hdul[0].header
        del header['*3']
        del header['*4']
        header['WCSAXES']= 2
        x_pos=header['CRVAL1']    
        y_pos=header['CRVAL2']    
        fits.writeto(new_fits_image, header=header, data=image_data, overwrite=True)
    return x_pos,y_pos,new_fits_image


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
    
def plot(fits_image1, fits_image2, name):
    # Load image data
    sigma=10
    vmin=1*sigma
    vmax=20*sigma
    x_pos,y_pos,fits_image1 = get_fitsimage(fits_image1,new_fits_image='image1_fix.fits')
    x_pos,y_pos,fits_image2= get_fitsimage(fits_image2,new_fits_image='image2_fix.fits')
    scale=scale_bar(z,pix_scale)
    fig = plt.figure(figsize=(18, 10))
    img1 = aplpy.FITSFigure(fits_image1, figure=fig,  subplot=[0.1,0.2,0.65,0.72])
    img1.recenter(x_pos,y_pos,radius=1.1)
    img1.axis_labels.set_font(size=18)
    img1.show_colorscale(vmin=vmin, vmax=vmax, cmap='inferno') 
    img1.axis_labels.set_xtext('RA (J2000)')
    img1.axis_labels.set_ytext('Dec (J2000)')
    img1.tick_labels.set_font(size=18)
    img1.ticks.set_length(10)
    img1.ticks.set_linewidth(2)

    img2 = aplpy.FITSFigure(fits_image2, figure=fig,  subplot=[0.55,0.1,0.35,0.7])
    img2.recenter(x_pos,y_pos,radius=0.76)
    img2.add_scalebar(scale)
    img2.scalebar.set_corner('bottom right')
    img2.scalebar.set_length(scale/3600)
    img2.scalebar.set_label(r'2 Mpc')
    img2.scalebar.set_color('black')
    img2.scalebar.set_font_size(size=20)
    img2.show_colorscale(vmin=vmin, vmax=vmax, cmap='inferno')  
    img2.axis_labels.set_xtext('RA (J2000)')
    img2.axis_labels.set_ytext('Dec (J2000)')
    img2.axis_labels.set_font(size=18)
    img2.tick_labels.set_font(size=18)
    img2.ticks.set_length(10)
    img2.ticks.set_linewidth(2)
    img2.add_beam()
    img2.beam.set_color('blue')
    img2.add_colorbar()
    img2.colorbar.set_axis_label_text(r'($\mu$Jy/Beam)')
    img2.colorbar.set_location('right')
    img2.colorbar.set_font(size=16)
    img2.colorbar.set_axis_label_font(size=26)
    img2.savefig('plots/A2631_plot_image.png', max_dpi=300) 
    # Set axis limits
  


if "__main__":
    z=0.27
    pix_scale=1.5
    name='A2631'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    plots='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'   
    fits_image1 = path+ 'mypipelinerun_ABELL2631_4-MFS-image.pbcor_full.fits'
    fits_image2 =path+'mypipelinerun_ABELL2631_4-MFS-image.pbcor.fits'
    # fits_image1 =path+ 'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.pbcor_full.fits'
    # fits_image2 =path+ 'mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.pbcor.fits'
    
    plot(fits_image1,fits_image2, name)
