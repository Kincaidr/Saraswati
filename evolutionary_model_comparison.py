from astropy.table import Table
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib.pyplot as plt
# Generate KDE curve for redshift
from matplotlib.ticker import FuncFormatter

def format_no_decimal(x, pos):
        if x == int(x):
            return str(int(x))
        else:
            return str(x)
# Apply to axes

TRECS_area=25
HSC_area=1.65
MIGHTEE_area=(1.04)**2

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
    Model_table = Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/catalogue_continuum_clustered_wrapped.fits')
    flux_model = Model_table['I1400']
    mask_flux = (flux_model < flux_max) & (flux_model > flux_min)
    mask_z = (Model_table['redshift'] > 0) & (Model_table['redshift'] < 4)
    mask_SFG = Model_table['logSFR'] != -100
    mask_AGN = Model_table['Lum1400'] != -100

    redshift_all = Model_table['redshift'][mask_flux & mask_z]
    redshift_SFR = Model_table['redshift'][mask_flux & mask_SFG & mask_z]
    redshift_AGN = Model_table['redshift'][mask_flux & mask_AGN & mask_z]
    return A2631_redshift_HSC,Zwcl2341_redshift_HSC, redshift_all, redshift_SFR, redshift_AGN

def redshift_table_Mightee(flux_min, flux_max):
    MIGHTEE_table = Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/MIGHTEE_catalog.fits')
    flux= MIGHTEE_table['S_INT14']*1e3  # Convert to mJy
    redshift= MIGHTEE_table['Redshift']
    mask_z = (redshift > 0) & (redshift < 4)
    mask_flux = (flux > flux_min) & (flux < flux_max)
    redshift_MIGHTEE = redshift[mask_flux & mask_z ]
    return redshift_MIGHTEE

flux_bins = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 1), (1, 10)]
fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharey=False)
bins=[15,15,10,10,10,8 ]
bins_curve=[8,8,8,8,8,6  ]
for i, (flux_min, flux_max) in enumerate(flux_bins):
    row, col = divmod(i, 3)
    ax = axes[row, col]
    A2631_redshift_HSC,Zwcl2341_redshift_HSC, redshift_all, redshift_SFR, redshift_AGN = redshift_table(flux_min, flux_max)
    MIGHTEE_Z=redshift_table_Mightee(flux_min, flux_max)
    print('Total N all',len(redshift_all)/TRECS_area)
    print('Total N A2631',len(A2631_redshift_HSC)/HSC_area)
    print('Total N Zwcl2341',len(Zwcl2341_redshift_HSC)/HSC_area)
    print('Total N MIGHTEE COSMOS',len(MIGHTEE_Z)/MIGHTEE_area)
    x1, y1 = np.histogram(A2631_redshift_HSC, bins=bins_curve[i], density=False, weights=np.ones_like(A2631_redshift_HSC) / HSC_area)
    x2, y2 = np.histogram(Zwcl2341_redshift_HSC, bins=bins_curve[i], density=False, weights=np.ones_like(Zwcl2341_redshift_HSC) / HSC_area)
    x3,y3 = np.histogram(MIGHTEE_Z, bins=bins_curve[i], density=False, weights=np.ones_like(MIGHTEE_Z) / MIGHTEE_area)
    bin_centers1 = 0.5 * (y1[1:] + y1[:-1])
    bin_centers2 = 0.5 * (y2[1:] + y2[:-1])
    bin_centers3 = 0.5 * (y3[1:] + y3[:-1])
    ax.plot(bin_centers1, x1, label='A2631', color='blue', linewidth=3)
    ax.plot(bin_centers2, x2, label='Zwcl2341', color='red', linewidth=3)
    ax.plot(bin_centers3, x3, label='MIGHTEE COSMOS', color='black', linewidth=3,linestyle='dashed')
    ax.hist(redshift_all, bins=bins[i], alpha=1, color='green',label='Total (AGN + SFG)',
            density=False, weights=np.ones_like(redshift_all) / TRECS_area, histtype='barstacked' ,rwidth=0.8)
    ax.hist(redshift_AGN, bins=bins[i], alpha=1, color='purple', label='AGN',
            density=False, weights=np.ones_like(redshift_AGN) / TRECS_area, histtype='barstacked', rwidth=0.8)
    
    ax.tick_params(axis='both', which='major', labelsize=12, length=5, width=1)
    ax.set_title(f'{flux_min}-{flux_max} mJy',fontsize=17)
    ax.text(0.98, 0.2, f'N = {len(redshift_all)}', transform=ax.transAxes,
        fontsize=15, ha='right', va='bottom')
     
    yticks = ax.get_yticks()
    ax.set_yticks(yticks.astype(int))

    if row ==1:
        ax.set_xlabel('Redshift',fontsize=17)
    if i == 0:
        ax.legend(fontsize=14)
    if col == 0:
        ax.set_ylabel('N (deg$^{-2}$)', fontsize=17)
    else:
        ax.set_ylabel('', fontsize=17)
plt.tight_layout()
plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/evolutionary_plots.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
plt.show()

