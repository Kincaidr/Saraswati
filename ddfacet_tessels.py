
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from astropy.wcs import WCS
from astropy.cosmology import Planck18 as cosmo
from astropy.wcs.utils import pixel_to_skycoord
from matplotlib.patches import Circle
from regions import Regions

def image_analyse(fits_image):

    hdu = fits.open(fits_image)[0]
    header = fits.getheader(fits_image)
    wcs = WCS(header)
    wcs=wcs[0,0,:,:]
    image_data = hdu.data

    if len(image_data.shape) > 2:
        image_data = image_data[0, 0,:, :] 
        
    sigma=8
    vmin=0.3*sigma
    vmax=10*sigma

    return(image_data,vmin, vmax, wcs)


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

def plot(fits_image1, name, regfile):
    # Load image data
    image_data1, vmin1, vmax1, wcs = image_analyse(fits_image1)
    fig, axes = plt.subplots(figsize=(12, 5 ), subplot_kw={'projection': wcs})  # Two side-by-side plots
    cax1 = axes.imshow(image_data1*1e6, cmap='grey_r', origin='lower', interpolation='none', vmin=vmin1, vmax=vmax1)
    axes.set_xlabel('RA (J2000)', size=18)
    axes.set_ylabel('Dec (J2000)', size=18)
    axes.tick_params(axis='both', which='major', labelsize=13, length=5, width=1)
    draw_ds9_lines(axes, regfile, color='black', linewidth=1)    
    plt.subplots_adjust(left=0.1, right=0.85, top=0.90, bottom=0.05) 
    plt.savefig('plots/'+name+'_tessel_image.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()


if "__main__":
    name='Zwcl2341'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    plots='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'   
    #fits_image = path+ 'image_DI_Clustered.DeeperDeconv.AP.app.restored.fits'
    fits_image =path+'image_DI_Clustered.DeeperDeconv.AP.app.restored.fits'
    # region_file= path+'image_DI_Clustered.tessel_facet_all.reg'
    region_file= path+name+'.tessel_facet_all.reg'
    plot( fits_image,name, region_file)
