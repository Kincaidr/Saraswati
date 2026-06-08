
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import numpy as np


name='A2631'
path=name+'/catalogs/'
cat1=name+'/catalogs/'+name+'_cut_srl.fits'
multi_comp_cat=path+'/multi_component_cat.fits'
table1 = Table.read(cat1)
table2 = Table.read(multi_comp_cat)

source_id_1 = table1['Source_id']
source_id_2 = table2['Source_id']
mask = np.isin(source_id_1,source_id_2)
table1=table1[~mask]
table1.write(path+name+'_multi_srl.fits',format='fits',overwrite=True)
