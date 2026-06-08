from astropy.table import Table
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
# Generate KDE curve for redshift

HSC_area=1.6
TRECS_area=25

def HSC_table(table):
    flux = table['Total_flux']*1e3  # Convert to mJy
    mask_flux_HSC = (flux > flux_min) & (flux < flux_max)
    mask_z_HSC = (table['photoz_best'] > 0) & (table['photoz_best'] < 4)
    redshift_HSC = table['photoz_best'][mask_z_HSC & mask_flux_HSC]
    return redshift_HSC

def redshift_table(flux_min, flux_max):
    HSC1_table = Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/A2631_HSC_dnnz_MeerKAT.fits')
    HSC2_table = Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/Zwcl2341_HSC_dnnz_MeerKAT.fits')
    A2631_redshift_HSC=HSC_table(HSC1_table)
    Zwcl2341_redshift_HSC=HSC_table(HSC2_table)
    Model_table = Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/catalogue_continuum_wrapped.fits')
    flux_model = Model_table['I1400']
    mask_flux = (flux_model < flux_max) & (flux_model > flux_min)
    mask_z = (Model_table['redshift'] > 0) & (Model_table['redshift'] < 4)
    mask_SFG = Model_table['logSFR'] != -100
    mask_AGN = Model_table['Lum1400'] != -100

    redshift_all = Model_table['redshift'][mask_flux & mask_z]
    redshift_SFR = Model_table['redshift'][mask_flux & mask_SFG & mask_z]
    redshift_AGN = Model_table['redshift'][mask_flux & mask_AGN & mask_z]
    return A2631_redshift_HSC,Zwcl2341_redshift_HSC, redshift_all, redshift_SFR, redshift_AGN


flux_bins = [(0.1, 0.2), (0.2, 0.4), (0.4, 1), (1, 2), (2, 20), (20, 40)]
fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharey=False)
bins=[35,35,35,15,15,8  ]
for i, (flux_min, flux_max) in enumerate(flux_bins):
    row, col = divmod(i, 3)
    ax = axes[row, col]
    
    A2631_redshift_HSC,Zwcl2341_redshift_HSC, redshift_all, redshift_SFR, redshift_AGN = redshift_table(flux_min, flux_max)
    x1, y1 = np.histogram(A2631_redshift_HSC, bins=20, density=False, weights=np.ones_like(A2631_redshift_HSC) / HSC_area)
    x2, y2 = np.histogram(Zwcl2341_redshift_HSC, bins=20, density=False, weights=np.ones_like(Zwcl2341_redshift_HSC) / HSC_area)
    bin_centers1 = 0.5 * (y1[1:] + y1[:-1])
    bin_centers2 = 0.5 * (y2[1:] + y2[:-1])
    ax.plot(bin_centers1, x1, label='A2631', color='blue', linewidth=3)
    ax.plot(bin_centers2, x2, label='Zwcl2341', color='red', linewidth=3)
    ax.hist(redshift_all, bins=bins[i], alpha=1, color='green',label='Total (AGN + SFG)',
            density=False, weights=np.ones_like(redshift_all) / TRECS_area, histtype='barstacked' ,rwidth=0.9)
    ax.hist(redshift_AGN, bins=bins[i], alpha=1, color='purple', label='AGN',
            density=False, weights=np.ones_like(redshift_AGN) / TRECS_area, histtype='barstacked', rwidth=0.9)

    ax.tick_params(axis='both', which='major', labelsize=12, length=5, width=1)
    if row ==1:
        ax.set_xlabel('Redshift',fontsize=16)
    ax.set_title(f'{flux_min}-{flux_max} mJy')
    if i == 0:
        ax.legend()
    
    if col == 0:
        ax.set_ylabel('N (deg$^{-2}$)',fontsize=16)
    else:
        ax.set_yticklabels([])
        ax.set_yticks([])

plt.tight_layout()
plt.show()

