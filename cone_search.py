from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

def cone_search(catalog,ra_center, dec_center, radius_arcmin):
    table = Table.read(catalog)
    catalog_coords = SkyCoord(ra=table['RA'], dec=table['DEC'], unit='deg')
    center_coord = SkyCoord(ra=ra_center, dec=dec_center, unit='deg')
    separation = center_coord.separation(catalog_coords)
    mask = separation < (radius_arcmin * u.arcmin)
    cone_table = table[mask]
    print(f"Number of objects in cone search: {len(cone_table)}")
    cone_table.write(path+name+'_multi_srl_new_cone.fits', format='fits', overwrite=True)
 
if __name__ == "__main__":
    name='Zwcl2341'
    name, cluster_centre='Zwcl2341',SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
    #name, cluster_centre='A2631',SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
    catalog=path+name+'_multi_srl_new.fits'
    ra_center=cluster_centre.ra 
    dec_center=cluster_centre.dec
    radius_arcmin=60
    cone_search(catalog,ra_center, dec_center, radius_arcmin)
    print('Cone search completed.')
    
