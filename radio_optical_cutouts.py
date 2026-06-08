import matplotlib.pyplot as plt
import requests
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import numpy as np
import os.path as op
import aplpy
from matplotlib.patches import Circle
import time

def get_fitsimage(fits_image):
    new_fits_image=name+'_deleted_axes.fits'
    if op.exists(new_fits_image): 
        print("You've already made this fits file")
        return new_fits_image
    with fits.open(fits_image) as hdul:
        image_data = hdul[0].data
        image_data=image_data[0,0,:,:]
        header=hdul[0].header
        del header['*3']
        del header['*4']
        header['WCSAXES']= 2    
        fits.writeto(new_fits_image, header=header, data=image_data, overwrite=True)
    return new_fits_image


def cutouts(table, fits_image,output_path):
    fits_image=get_fitsimage(fits_image)
    width,height=180,180

    for i, (ra,dec,ra1,dec1,source_id1,source_id2) in enumerate(zip(table['RA'],table['DEC'], table['RA2'],table['DEC2'], table['Source_id'], table['Source_id_2'])):
        print(f"Processing entry #{i}")
        coord1 = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        coord2 = SkyCoord(ra=ra1, dec=dec1, unit=(u.deg, u.deg), frame='icrs') 
        center = SkyCoord(ra=(coord1.ra + coord2.ra) / 2,dec=(coord1.dec + coord2.dec) / 2)  
        image_name=output_path+ f'rgb_cutout_{source_id1}_{source_id2}.jpg'
        optical_name=output_path+f'optical_cutout_{source_id1}_{source_id2}.fits'
        optical_image = get_decals(optical_name, center, width, height)
        time.sleep(2)
        aplpy.make_rgb_image(optical_image, image_name,stretch_b='linear',stretch_g='linear',stretch_r='linear',pmax_b=99,pmax_g=99,pmax_r=99)
        fig = plt.figure(figsize=(5, 5))
        img = aplpy.FITSFigure(image_name, figure=fig)
        img.show_rgb()
        levels = np.array([2,4,8,16,32])*sigma
        img.show_contour(fits_image, colors='lightblue', linewidths=0.6, alpha=1.0, levels=levels, smooth=1)
        img.show_circles(ra, dec, radius=0.0003, edgecolor='lime', lw=1)
        img.show_circles(ra1, dec1, radius=0.0003, edgecolor='lime', lw=1)
        #img.show_circles(ra_gaia, dec_gaia, radius=0.0003, edgecolor='blue', lw=1)
        img.axis_labels.set_font(size=15)
        img.add_scalebar(20/3600)
        img.scalebar.set_length(20/3600)
        img.scalebar.set_label(r'$20^{\prime\prime}$')
        img.scalebar.set_color('yellow')
        img.scalebar.set_font_size(size=15)

        #img.recenter(center.ra,center.dec, width =0.04,height = 0.04)
        img.savefig(image_name, dpi=300)
        plt.close(fig)  # ✅ Frees memory
        del img
        del fig
        print("Saved cutout image:", image_name)


def get_decals(file_name, pos, width,height, pixmax=3000, justone=None):    
    pixscale = 0.262 
    if op.exists(file_name): 
        print("You've already made this file")
        return file_name
    fname = 'cutout.fits?ra={}&dec={}&layer=hsc-dr3&pixscale={}&width={}&height={}'.format(pos.ra.deg, pos.dec.deg, pixscale,width, height)
    
    url = 'https://www.legacysurvey.org/viewer/{}'.format(fname)
    #url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=150'    
    #url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=150'    
    response = requests.get(url)
    print("url retreived:", url)
    if response.status_code == 200:
        with open(file_name, "wb") as f:
            f.write(response.content)
        print("Image saved:", file_name)
    return file_name
            

if __name__ == "__main__":
    name='Zwcl2341'
    sigma=11e-6
    catalog_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/multi_component_cat.fits'
    fits_image='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/mypipelinerun_ABELL2631_4-MFS-image.pbcor.fits'#mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.pbcor.fits'
    output_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/cutouts/'
    table=Table.read(catalog_path)
    chunk_size=10

    for start in range(0,len(table),chunk_size):
        end = min(start + chunk_size, len(table))
        print(f"Processing entries {start} to {end-1}")
        cutouts(table[start:end], fits_image, output_path)
        print('start',start)
        #Wait 10 seconds after each chunk (except the last one)
        if end < len(table):
            print("Waiting 20 seconds before next chunk...")
            time.sleep(10)


    