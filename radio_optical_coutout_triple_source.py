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


def cutouts(fits_image,output_path,catalog):
    fits_image=get_fitsimage(fits_image)
    width,height=250,250
    # source_ids1=[1152,1153]
    # source_ids2=[1152,1154]
    #source_ids1=[875,876]
    source_ids1 =[345,346]
    source_ids2=[346,347]
    # source_ids1 =[633,635]
    # source_ids2=[634,635]
    #source_ids1=[1161,1162]
    #source_ids2=[1161,1163]
    # #source_ids1=[1974,1975]
    # #source_ids2=[1974,1976]
    # source_ids1=[1765,1767]
    # source_ids2=[1765,1766]
    # source_ids1=[1618,1619]
    # source_ids2=[1618,1620]
    # source_ids1=[918,920]
    # source_ids2=[919,920]

    table=Table.read(catalog)
    source_id=table['Source_id']  
    source_id2=table['Source_id_2']  
    mask1=(source_id==source_ids1[0])&(source_id2==source_ids1[1])
    mask2=(source_id==source_ids2[0])&(source_id2==source_ids2[1])
    ra,dec=table['RA'][mask1],table['DEC'][mask1]
    ra1,dec1=table['RA2'][mask1],table['DEC2'][mask1]
    ra2,dec2=table['RA'][mask2],table['DEC'][mask2]
    ra3,dec3=table['RA2'][mask2],table['DEC2'][mask2]
    breakpoint()
    coord1 = SkyCoord(ra=ra1, dec=dec1, unit=(u.deg, u.deg), frame='icrs')
    image_name=output_path+ f'rgb_cutout_triple_source_{source_ids1[0]}_{source_ids1[1]}_{source_ids2[1]}.jpg'
    optical_name=output_path+f'optical_cutout_{source_ids1[0]}_{source_ids1[1]}_{source_ids2[1]}.fits'
    optical_image = get_decals(optical_name, coord1, width, height)
    aplpy.make_rgb_image(optical_image, image_name,stretch_b='linear',stretch_g='linear',stretch_r='linear',pmax_b=99,pmax_g=99,pmax_r=99)
    fig = plt.figure(figsize=(5, 5))
    img = aplpy.FITSFigure(image_name, figure=fig)
    img.show_rgb()
    levels = np.array([2,4,8,16,32])*sigma
    img.show_contour(fits_image, colors='lightblue', linewidths=0.6, alpha=1.0, levels=levels, smooth=1)
    img.show_circles(ra, dec, radius=0.0003, edgecolor='lime', lw=1)
    img.show_circles(ra1, dec1, radius=0.0003, edgecolor='lime', lw=1)
    img.show_circles(ra2, dec2, radius=0.0003, edgecolor='lime', lw=1)
    img.show_circles(ra3, dec3, radius=0.0003, edgecolor='lime', lw=1)
    img.axis_labels.set_font(size=15)
    img.add_scalebar(20/3600)
    img.scalebar.set_length(20/3600)
    img.scalebar.set_label(r'$20^{\prime\prime}$')
    img.scalebar.set_color('yellow')
    img.scalebar.set_font_size(size=15)
    img.savefig(image_name,dpi=300)
    plt.close(fig)  
    print("Saved cutout image:", image_name)


def get_decals(file_name, pos, width,height, pixmax=3000, justone=None):
    # Code from SoFiA-image-pipeline
    # Get DECaLS false color image and fits (for the WCS). Example URL for this script provided by John Wu.
    pixscale = 0.262 
    pixscale = 0.3 
    if op.exists(file_name): 
        print("You've already made this file")
        return file_name
    fname = 'cutout.fits?ra={}&dec={}&layer=hsc-dr3&pixscale={}&width={}&height={}'.format(pos.dec.deg, pos.ra.deg, pixscale,width, height)
    
    url = 'https://www.legacysurvey.org/viewer/{}'.format(fname)
    #url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=150'    
    #url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=150'    
    response = requests.get(url)
    print("url retreived:", url)

    with open(file_name, "wb") as f:
        f.write(response.content)
    print("Image saved:", file_name)
    return file_name
            

if __name__ == "__main__":
    name='A2631'
    sigma=10e-6
    catalog='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/multi_component_cat.fits'
    fits_image='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/mypipelinerun_ABELL2631_4-MFS-image.pbcor.fits'
    #fits_image='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.pbcor.fits'
    output_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/cutouts/'
    cutouts(fits_image,output_path,catalog)




    