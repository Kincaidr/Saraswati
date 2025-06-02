
import matplotlib.pyplot as plt
import requests
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import numpy as np
import os.path as op


def get_decals(file_name, pos, width,height, pixscale=0.262):
    # Code from SoFiA-image-pipeline
    # Get DECaLS false color image and fits (for the WCS). Example URL for this script provided by John Wu.

    if op.exists(file_name): 
        print("You've already made this file")
        return file_name
    fname = 'cutout.fits?ra={}&dec={}&layer=hsc-dr3&pixscale={}&width={}&height={}'.format(pos.ra.deg, pos.dec.deg, pixscale,width, height)
    
    #url = 'https://www.legacysurvey.org/viewer/{}'.format(fname)
    url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=10000'    
    #url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra='+str(pos.ra.value)+'&dec='+str(pos.dec.value)+'&layer=hsc-dr3&&pixscale=0.26&size=150'    
    response = requests.get(url)
    print("url retreived:", url)
    
    with open(file_name, "wb") as f:
            f.write(response.content)
    # if response.status_code == 200:
    #     with open(file_name, "wb") as f:
    #         f.write(response.content)
    #     print("Image saved:", file_name)
    return file_name


if __name__ == "__main__":
    name='A2631' 
    ra,dec=-5.5808,0.2766
    width,height = 20000,20000
    pos=SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    get_decals(name+'_HSC_full.fits', pos, width,height)