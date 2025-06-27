import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy import stats
from astropy import units as u
from astropy.coordinates import SkyCoord

def scale_flux(flux, counts_freq=1283, data_freq=1400, spectral_index=-0.7):
    return flux * (counts_freq / data_freq) ** spectral_index

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
        peak1 = cat1['Peak_flux']*1e3
        peak2 = cat2['Peak_flux']*1e3
        flux2= scale_flux(flux2, counts_freq=1283, data_freq=1400, spectral_index=-0.7)
        snr1=cat1['Isl_rms']*1e3
        snr2=cat2['Isl_rms']*1e3
        mask= (flux1/snr1 > 20) & (flux2/snr2 > 20) 
        #mask2= (flux1/peak1 > 0.8) & (flux1/peak1 < 1.2) & (flux2/peak2 > 0.8)  & (flux2/peak2 < 1.2) 

        table1=cat1[mask]
        table2=cat2[mask]
        flux1 = table1['Total_flux']*1e3
        flux2 = table2['Total_flux']*1e3
        flux1_err = table1['E_Total_flux']*1e3
        flux2_err = table2['E_Total_flux']*1e3
        print('Total number of sources',len(flux1))
        print('Mean flux ratio',np.mean(flux1/flux2))
        print('Standard deviation ratio',np.std(flux1/flux2))
        xx = np.linspace(1,1000,1024)
        fig = plt.figure(figsize=(15, 12))
        plt.plot(xx,xx, "k--",linewidth=4, alpha=0.7)   
        #plt.scatter(flux1, flux2, s=50, c='orange', alpha=0.5, edgecolors='black', linewidth=1)
        slope, intercept,r_value, p_value, std_err = stats.linregress(flux1, flux2)
        print('r_value',r_value)
        plt.plot(xx,linear_function(xx, slope, intercept),linewidth=4, alpha=0.5)
        plt.errorbar(flux1, flux2, xerr=flux1_err, yerr=flux2_err, fmt='o', color="#D19C28", alpha=0.5, ecolor='black', elinewidth=1, capsize=2)
        plt.ylabel(r'MeerKAT $S_{T}$ [mJy]',fontsize=26)
        plt.xlabel(r'FIRST $S_{T}$ [mJy]',fontsize=26)
        plt.tick_params(axis='both', which='major', direction='in', length=10, width=2, labelsize=22)
        plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=22)
        plt.xscale('log')
        plt.yscale('log')
        axtop = plt.axes([0.125, 0.85, 0.75, 0.15])
        axtop.hist(flux1, bins=20, range=(0,50),color='blue', alpha=0.7, orientation='vertical',histtype='step',linewidth=4)
        axtop.axvline(np.median(flux1),linestyle='dashed',color='black',linewidth=2)
        axright = plt.axes([0.87, 0.11, 0.13, 0.75])
        axright.hist(flux2 , bins=20,range=(0,50), color='blue', alpha=0.7, orientation='horizontal',histtype='step',linewidth=4)
        axright.axhline(np.median(flux1  ),linestyle='dashed',color='black',linewidth=2)
        axtop.set_xticks([])
        axtop.set_yticks([])
        axright.set_xticks([])
        axright.set_yticks([])
        axright.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
        axtop.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
        plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/FIRST_flux_scale.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

if __name__ == "__main__":
    name='Zwcl2341'  # Change this to the desired cluster name
    path=f'/home/kincaid/Desktop/Saraswati_codes/'
    cat1=f'{path}/MeerKAT_combined_new.fits'
    cat2=f'{path}/FIRST_combined.fits'
    table1 = Table.read(cat1)
    table2 = Table.read(cat2)
    cat1_matched,cat2_matched=crossmatch_catalogs(table1, table2, radius=5)
    plot(cat1_matched,cat2_matched)
