import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy import stats
from astropy import units as u
from astropy.coordinates import SkyCoord

def linear_function(x, m, b):
        return m * x + b

def crossmatch_catalogs(cat1, cat2, radius=5):
    coords1 = SkyCoord(cat1['RA'], cat1['DEC'], unit=(u.deg, u.deg))
    breakpoint()
    coords2 = SkyCoord(cat2['RA'], cat2['DEC'],unit=(u.deg, u.deg))
    idx, d2d, _ = coords1.match_to_catalog_sky(coords2)
    matched = d2d.arcsec < radius  # Apply matching radius condition
    return cat1[matched], cat2[idx[matched]]

def plot(cat1,cat2):
        flux1 = cat1['Total_flux']*1e3
        flux2 = cat2['Total_flux']*1e3
        snr1=cat1['Isl_rms']*1e3
        snr2=cat2['Isl_rms']*1e3
        mask= (flux1/snr1 > 20) & (flux2/snr2 > 20)
        table1=cat1[mask]
        table2=cat2[mask]
        flux1 = table1['Total_flux']*1e3
        flux2 = table2['Total_flux']*1e3
        print('Total number of sources',len(flux1))
        print('Mean flux ratio DDFacet wsclean',np.mean(flux1/flux2))
        xx = np.linspace(0,1000,1024)
        fig = plt.figure(figsize=(15, 12))
        plt.plot(xx,xx, "k--",linewidth=4, alpha=0.5)   
        plt.scatter(flux1, flux2, s=50, c='orange', alpha=0.5, edgecolors='black', linewidth=1)
        plt.ylabel(r'DDFacet PB corr $S_{T}$ [mJy]',fontsize=26)
        plt.xlabel(r'WSClean katbeam PB corr $S_{T}$ [mJy]',fontsize=26)
        plt.tick_params(axis='both', which='major', direction='in', length=10, width=2, labelsize=22)
        plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=22)
        plt.xscale('log')
        plt.yscale('log')
        axtop = plt.axes([0.125, 0.85, 0.75, 0.15])
        axtop.hist(flux1, bins=20, range=(0,50),color='orange', alpha=0.7, orientation='vertical',histtype='step',linewidth=4)
        axtop.axvline(np.median(flux1),linestyle='dashed',color='black',linewidth=2)
        axright = plt.axes([0.87, 0.11, 0.13, 0.75])
        axright.hist(flux2 , bins=20,range=(0,50), color='orange', alpha=0.7, orientation='horizontal',histtype='step',linewidth=4)
        axright.axhline(np.median(flux1  ),linestyle='dashed',color='black',linewidth=2)
        axtop.set_xticks([])
        axtop.set_yticks([])
        axright.set_xticks([])
        axright.set_yticks([])
        axright.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
        axtop.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
        plt.savefig('FIRST_flux_scale.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

if __name__ == "__main__":
    cat1='A2631/catalogs/A2631_srl.fits'
    cat2='A2631/catalogs/A2631_katbeam_srl.fits'
    table1 = Table.read(cat1)
    table2 = Table.read(cat2)
    cat1_matched,cat2_matched=crossmatch_catalogs(table1, table2, radius=5)
    plot(cat1_matched,cat2_matched)
