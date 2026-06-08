import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy import stats
from astropy import units as u
from astropy.coordinates import SkyCoord

def linear_function(x, m, b):
        return m * x + b

def plot(cat1,freq1,freq2):
        table=Table.read(cat1)
        flux1 = table['Total_flux']*1e3
        flux2 = table['Fint']
        peak1=table['Peak_flux']*1e3
        peak2=table['Fpeak']  
        rms=table['Isl_rms']*1e3
        flux1_err=table['E_Total_flux']*1e3
        RA1,DEC1=table['RA_1'],table['DEC_1']
        RA2,DEC2=table['RA_2'],table['DEC_2']
        coord1 = SkyCoord(ra=RA1, dec=DEC1, unit=(u.deg, u.deg), frame='icrs')  
        coord2 = SkyCoord(ra=RA2, dec=DEC2, unit=(u.deg, u.deg), frame='icrs')
        sep = coord1.separation(coord2).arcsec  
        mask= (flux1/rms > 10) & (flux1/peak1 > 0.8) & (flux2/peak2 > 0.8)  & (flux1/peak1 < 1.2) & (flux2/peak2 < 1.2) 
        table=table[mask]
        flux1 = table['Total_flux']*1e3
        flux2 = table['Fint']
        peak1=table['Peak_flux']*1e3
        peak2=table['Fpeak']  
        flux1_err=table['E_Total_flux']*1e3
        print('Total number of sources',len(flux1))
        flux2_err = np.ones(len(flux2)) * 0.1
        flux1_corr= flux1* (freq2 / freq1)**-0.7
        flux1_corr_err= flux1_err * (freq2 / freq1)**-0.7
        flux2_corr = flux2
        flux2_corr_err= flux2_err 

        flux_ratio=flux2_corr/flux1_corr
        print('Mean flux ratio',np.mean(flux_ratio))
        print('Standard deviation ratio',np.std(flux_ratio))
        correction_factor=1+(1-np.mean(flux_ratio))
        print('Correction factor',correction_factor)
        flux2_corr=flux2_corr*(correction_factor)
        flux_ratio=flux2_corr/flux1_corr
        print('New Mean flux ratio',np.mean(flux_ratio))
        print('New Standard deviation ratio',np.std(flux_ratio))
        xx = np.linspace(0,1000,1024)
        fig = plt.figure(figsize=(15, 12))
        plt.plot(xx,xx, "k--",linewidth=4, alpha=0.5)   
        slope, intercept,r_value, p_value, std_err = stats.linregress(np.log10(flux1_corr), np.log10(flux2_corr))
        print('r_value',r_value)
        plt.plot(xx,linear_function(xx, slope, intercept),linewidth=4, alpha=0.5)
        plt.errorbar(flux1_corr,flux2_corr ,
                        yerr=abs(flux1_corr_err),
                        xerr=abs(flux2_corr_err),
                        capsize=7, linewidth=1,
                        color="k",markersize=2,marker='o', alpha=0.7,linestyle='None')
        plt.xlabel(r'MeerKAT $S_{T}$ [mJy]',fontsize=26)
        plt.ylabel(r'FIRST $S_{T}$ [mJy]',fontsize=26)
        plt.tick_params(axis='both', which='major', direction='in', length=10, width=2, labelsize=22)
        plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1.5, labelsize=22)
        plt.xscale('log')
        plt.yscale('log')
        axtop = plt.axes([0.125, 0.85, 0.75, 0.15])
        axtop.hist(flux1_corr, bins=20, range=(0,50),color='orange', alpha=0.7, orientation='vertical',histtype='step',linewidth=4)
        axtop.axvline(np.median(flux1_corr),linestyle='dashed',color='black',linewidth=2)
        axright = plt.axes([0.87, 0.11, 0.13, 0.75])
        axright.hist(flux2_corr , bins=20,range=(0,50), color='orange', alpha=0.7, orientation='horizontal',histtype='step',linewidth=4)
        axright.axhline(np.median(flux1_corr  ),linestyle='dashed',color='black',linewidth=2)
        axtop.set_xticks([])
        axtop.set_yticks([])
        axright.set_xticks([])
        axright.set_yticks([])
        axright.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
        axtop.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
        plt.savefig('FIRST_flux_scale.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

if __name__ == "__main__":
    freq1 = 1.283
    freq2 = 1.4  
    cat1='catalogs/A2631_Zwcl_combined_ddfacet_FIRST.fits'
    plot(cat1,freq1,freq2)
