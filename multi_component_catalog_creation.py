
import numpy as np
from astropy.table import Table, Column
from astropy import units as u
from astropy.coordinates import SkyCoord

name='A2631'
catalog_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
catalog=catalog_path+'multi_component_cat.fits'
table=Table.read(catalog)

RA=table['RA']
DEC=table['DEC']
RA2=table['RA2']
DEC2=table['DEC2']
Source_id=table['Source_id']
Source_id2=table['Source_id_2']

RA_all=Column(list(RA)+list(RA2))
DEC_all=Column(list(DEC)+list(DEC2))
Source_id_all=Column(list(Source_id)+list(Source_id2))
multi_component=Column(list(Source_id2)+list(Source_id))
 
catalog = Table({'RA': RA_all, 'DEC': DEC_all,'Source_id':Source_id_all, 'Multi_component':multi_component})
catalog.write(catalog_path+'multi_component_cat_real.fits', format='fits', overwrite=True)

