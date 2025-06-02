from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def script(cat1, name):

    table=Table.read(cat1)
    table.rename_column('Isl_rms', 'RMS')
    table.rename_column('Source_id', 'ID')
    RA=table['RA']+360
    DEC=table['DEC']
    table['RA']=table['RA']+360
    table['Total_flux']=table['Total_flux']*1e6
    table['E_Total_flux']=table['E_Total_flux']*1e6
    table['Peak_flux']=table['Peak_flux']*1e6
    table['E_Peak_flux']=table['E_Peak_flux']*1e6
    table['RMS']=table['RMS']*1e6
    table['Maj']=table['Maj']*3600
    table['Min']=table['Min']*3600

    cols = table.colnames
    coord = SkyCoord(ra=RA, dec=DEC, unit=(u.deg, u.deg), frame='icrs') 
    ra_str = coord.ra.to_string(unit=u.hourangle, sep='', pad=True, precision=2)    
    dec_str = coord.dec.to_string(unit=u.deg, sep='', alwayssign=True, pad=True, precision=1) 

    table['Name'] = [f"{name} J{ra}{dec}" for ra, dec in zip(ra_str, dec_str)]

    print(table['Name'])

    new_order =['Name'] + [col for col in cols if col != ['ID','Name']]
    table = table[new_order]
    new_cat=cat1=path+name+'_final_srl.fits'
    table.write(new_cat, overwrite=True)

    table['RA']=[f"{val:.6f}" for val in table['RA']]
    table['E_RA']=[f"{val:.5f}" for val in table['E_RA']]
    table['DEC']=[f"{val:.6f}" for val in table['DEC']]
    table['E_DEC']=[f"{val:.5f}" for val in table['E_DEC']]
    table['Total_flux']=[f"{val:.0f}" for val in table['Total_flux']]
    table['E_Total_flux']=[f"{val:.2f}" for val in table['E_Total_flux']]
    table['Peak_flux']=[f"{val:.0f}" for val in table['Peak_flux']]
    table['E_Peak_flux']=[f"{val:.2f}" for val in table['E_Peak_flux']]
    table['RMS']=[f"{val:.2f}" for val in table['RMS']]
    table['Maj']=[f"{val:.2f}" for val in table['Maj']]
    table['Min']=[f"{val:.2f}" for val in table['Min']]
    print(table['ID','Name','RA','E_RA','DEC','E_DEC','Maj','Min','Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','RMS','S_Code'])
    breakpoint()

if __name__ == "__main__":

    name='Zwcl2341'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
    cat1=path+name+'_srl_flux_corr.fits'
    script(cat1, name)    
    