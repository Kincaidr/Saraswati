
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord

def multi_component(catalog_path,catalog):
    
    table=Table.read(catalog)
    S_code=table['S_Code']
    mask=(S_code=='M') | (S_code=='C') 
    table_new=table[mask]
    S_code=table_new['S_Code']
    RA=table_new['RA']
    DEC=table_new['DEC']
    flux=table_new['Total_flux']
    Source_id=table_new['Source_id']
    breakpoint()
    x=len(table_new['S_Code'])
    print('Number of sources with S_Code =M,C',x)
    print('Percentage of sources with S_Code =M,,C',(x/len(table['S_Code']))*100)
    
    pos = SkyCoord(RA, DEC, frame='icrs', unit=(u.deg, u.deg))
    flux_ratio_matrix = flux[:, None] / flux[None, :]  
    sum_flux_matrix = flux[:, None] + flux[None, :] 
    sep_matrix= pos[:, None].separation( pos[None, :]).arcsec
    mask = ((flux_ratio_matrix >= 1/6) & (flux_ratio_matrix <= 6)).astype(int)
    sum_flux_matrix=sum_flux_matrix*mask
    sep_matrix=sep_matrix*mask
    theta_matrix = 100*np.sqrt(sum_flux_matrix*1e3/20) 
    row_indices, col_indices=np.where(sep_matrix  < theta_matrix)
    mask = row_indices != col_indices
    filtered_rows = row_indices[mask]
    filtered_cols = col_indices[mask]
    true_pairs = np.column_stack((filtered_rows, filtered_cols))
    positions=pos[np.unique(np.sort(true_pairs,axis=1),axis=0)]

    all_ra1=[]
    all_dec1=[]
    all_ra2=[]
    all_dec2=[]
    for pos in positions:
        ra1=pos[0].ra.value
        dec1=pos[0].dec.value
        ra2=pos[1].ra.value
        dec2=pos[1].dec.value
        all_ra1.append(ra1)
        all_dec1.append(dec1)
        all_ra2.append(ra2)
        all_dec2.append(dec2)

    pos = SkyCoord(RA, DEC, frame='icrs', unit=(u.deg, u.deg))
    pair1_coords = SkyCoord(ra=np.array(all_ra1)*u.deg, dec=np.array(all_dec1)*u.deg)
    pair2_coords = SkyCoord(ra=np.array(all_ra2)*u.deg, dec=np.array(all_dec2)*u.deg)
    idx1, sep1, _ = pair1_coords.match_to_catalog_sky(pos)
    idx2, sep2, _ = pair2_coords.match_to_catalog_sky(pos)

    source_id_1 = Source_id[idx1]
    source_id_2 = Source_id[idx2]
    S_Code_1 =  S_code[idx1]
    S_Code_2 =  S_code[idx2]
    catalog = Table({'RA': np.array(all_ra1), 'DEC': np.array(all_dec1),'Source_id':np.array(source_id_1), 'S_Code':np.array(S_Code_1), 'RA2': np.array(all_ra2), 'DEC2': np.array(all_dec2),'Source_id_2':np.array(source_id_2), 'S_Code_2':np.array(S_Code_2)})
    catalog.write(catalog_path+'multi_component_cat.fits', format='fits', overwrite=True)
    print("Number of multi-component sources:", len(catalog))
    return mask


if "__main__":
    name='A2631'
    catalog_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
    catalog=catalog_path+name+'_resolution_corr_srl.fits'
    num=multi_component(catalog_path,catalog)
