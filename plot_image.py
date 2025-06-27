
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from astropy.cosmology import Planck18 as cosmo
from astropy.wcs.utils import pixel_to_skycoord
from matplotlib.patches import Circle
from regions import Regions


def draw_ds9_lines(ax, regfile, color='black', linewidth=2):
    """
    Reads a DS9 .reg file with line(x1,y1,x2,y2) regions in pixel coordinates
    and plots them as solid lines on the given matplotlib Axes.
    """
    with open(regfile, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if line.startswith('line(') and ')' in line:
            coords_str = line.split('line(')[1].split(')')[0]
            try:
                x1, y1, x2, y2 = map(float, coords_str.split(','))
                ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle='-')
            except ValueError:
                print(f"Skipping invalid line: {line}")

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
        image_data = image_data[0, 0,:, :] 
        
    sigma=10
    vmin=1*sigma
    vmax=20*sigma
    return(image_data,vmin, vmax, wcs)

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
    image_data1, vmin1, vmax1, wcs = image_analyse(fits_image1)
    image_data2, vmin2, vmax2, _= image_analyse(fits_image2)
    scale=scale_bar(z,pix_scale)
    radius_arcsec = 0.8*3600
    radius_pix=  radius_arcsec/pix_scale# radius in degrees for the circle
    xpos, ypos = 1500, 1500
    fig, axes = plt.subplots(1, 2, figsize=(16, 7 ), subplot_kw={'projection': wcs})  # Two side-by-side plots
    cax1 = axes[0].imshow(image_data1 * 1e6, cmap='inferno', origin='lower', interpolation='none', vmin=vmin1, vmax=vmax1)
    axes[0].set_xlabel('RA (J2000)', size=16)
    axes[0].set_ylabel('Dec (J2000)', size=16)
    axes[0].set_xlim(500, 6000)   # "zoomed out"
    axes[0].set_ylim(500, 6000)
    #draw_ds9_lines(axes[0], regfile, color='white', linewidth=0.2)

    # axes[0].plot([xpos, xpos + scale], [ypos, ypos], color='yellow', lw=3)
    # axes[0].text(xpos + scale / 2, ypos + 50, '2 Mpc', color='white', fontsize=12, ha='center')
    # Plot second image
    cax2 = axes[1].imshow(image_data2 * 1e6, cmap='inferno', origin='lower', interpolation='none', vmin=vmin2, vmax=vmax2)
    axes[1].set_xlabel('RA (J2000)', size=16)
    axes[1].set_ylabel('Dec (J2000)', size=16)
    #axes[1].set_xlim(1200, 4800)   
    #axes[1].set_ylim(1200, 4800)
    axes[1].set_xlim(1200, 4800)
    axes[1].set_ylim(1100, 4900)
    axes[1].plot([xpos + 50, xpos + 50 + scale], [ypos + 50, ypos + 50], color='black', lw=4)
    axes[1].text(xpos + 50 + scale / 2, ypos + 120, '2 Mpc', color='black', fontsize=14, ha='center')

#     # Set axis limits
    for ax in axes:
        ax.tick_params(axis='both', which='major', labelsize=13, length=5, width=1)

    # if regfile:
    #     regions = Regions.read(regfile, format='ds9')
    #     for region in regions:
    #         patch = region.to_mpl(transform=axes[0].get_transform('world'))
    #         patch.set_edgecolor('black')     # <-- Change color here
    #         patch.set_linewidth(2.5)
    #         patch.set_facecolor('none')       # No fill
    #         patch.set_alpha(1.0)  
    #         axes[0].add_patch(patch)
    #         axes[1].add_patch(patch) 
    # # Single colorbar
    cbar_ax = fig.add_axes([0.88, 0.12, 0.02, 0.7])  # [left, bottom, width, height]
    colorbar = fig.colorbar(cax1, cax=cbar_ax, orientation='vertical')  # Attach colorbar to first plot
    colorbar.set_label(r'($\mu$Jy/Beam)', fontsize=16, color='black')
    plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to fit colorbar
    plt.subplots_adjust(left=0.1, right=0.85, top=0.90, bottom=0.05) 
    plt.savefig(plots+name+'_image.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()

if "__main__":
    z=0.27
    pix_scale=1.5
    name='A2631'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    plots='/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/'
    fits_image1 = path+ 'image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits'
    fits_image2 =path+ 'image_DI_Clustered.DeeperDeconv.AP4.int.restored.pbcor.fits'
    plot(fits_image1, fits_image2, name)
