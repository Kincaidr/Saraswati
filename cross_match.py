
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.table import Table

def cross_match(tab1,tab2, radius):
    table1=Table.read(tab1)
    table2=Table.read(tab2)
    ra1,ra2 = table1['RA'], table2['RA']
    dec1,dec2 = table1['DEC'], table2['DEC']
    cat1 = SkyCoord(ra=ra1, dec=dec1)
    cat2 = SkyCoord(ra=ra2, dec=dec2)
    idx, d2d, d3d = match_coordinates_sky(cat1, cat2)
    max_sep = radius * u.arcsec
    mask = d2d < max_sep
    table1 = table1[~mask]
    table2 = table2[idx[mask]]
    table1.write(f'{name}_cut_srl_filter.fits', overwrite=True)
    table2.write(f'{name}_filtered_table2.fits', overwrite=True)    
    return table1, table2


if __name__ == "__main__":
    name='Zwcl2341'  # Change this to the desired cluster name
    path=f'/home/kincaid/Desktop/Saraswati_codes/{name}/catalogs/'
    tab1=path+f'{name}_cut_srl.fits'
    tab2=path+f'{name}_res_srl.fits'
    cross_match(tab1,tab2, radius=5)
