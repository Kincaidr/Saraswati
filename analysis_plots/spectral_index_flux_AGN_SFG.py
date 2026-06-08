from astropy.table import Table
from scipy.stats import linregress
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u


def radio_luminosity(flux, z, alpha):
    flux_jy = flux * u.mJy  # observed flux density in Janskys
    flux = flux_jy.to(u.W / u.m**2 / u.Hz)
    L=4*np.pi*(cosmo.luminosity_distance(z).to(u.m))**2/((1+z)**(1+alpha))*((1.4/1.28)**alpha)*flux
    return L.value

def find_spectral_index(real_cat):
    cat = Table.read(real_cat)
    freqs = [998, 1283, 1569]  # Using the three frequencies
    spectral_index = []
    fluxes=[]
    source_ids=[]
    peaks=[]
    for i in range(len(cat['Total_flux'])):
        all_fluxes = [cat['Total_flux_1'][i], cat['Total_flux_2'][i], cat['Total_flux'][i]]
        all_freqs = [freqs[0], freqs[1], freqs[2]]
        Result = linregress(np.log10(all_freqs), np.log10(all_fluxes))
        source_ids.append(cat['Source_id'][i])
        spectral_index.append(Result[0])
        fluxes.append(cat['Total_flux_2'][i]*1e3)
        peaks.append(cat['Peak_flux_2'][i]*1e3)
    return(peaks,source_ids,fluxes,spectral_index)

def median_SI(fluxes, SI,nbins):
    #Range_x = 10**np.linspace(start=np.log10(fluxes.min()), stop=np.log10(fluxes.max()), num=nbins+1)
    Range_x = np.percentile(fluxes, np.linspace(0, 100, nbins + 1))
    bin_centers=[]
    median_sizes=[]
    errors=[]
    difflr=np.diff(Range_x)
    for i in range(len(Range_x)-1):
        mask = (fluxes >= Range_x[i]) & (fluxes < Range_x[i+1])
        bin_centers.append(np.sqrt(Range_x[i] * Range_x[i+1]))
        median_sizes.append(np.median(SI[mask]))
        error = 1.253 * np.std(SI[mask]) / np.sqrt(len(SI[mask]))
        errors.append(error)
    return(difflr,bin_centers, median_sizes, errors)

def plot(peakes,ids,fluxes,SI,AGN_SFR_NED_cat):
    cat2=Table.read(AGN_SFR_NED_cat)
    total_SI=np.concatenate(SI)
    total_fluxes=np.concatenate(fluxes)
    total_ids=np.concatenate(ids)
    total_peaks=np.concatenate(peakes)

    mask=total_fluxes < 200
    total_SI=np.array(total_SI[mask])
    total_fluxes=np.array(total_fluxes[mask])
    total_ids=np.array(total_ids[mask])
    total_peaks = np.array(total_peaks[mask])  
    common_ids,  idx_cat1, idx_cat2 = np.intersect1d(total_ids, cat2['Source_id'], return_indices=True)
    total_fluxes=total_fluxes[idx_cat1] 
    total_SI=total_SI[idx_cat1]
    total_peaks=total_peaks[idx_cat1]
    #z=cat2['Redshift (z)'][idx_cat2]
    z=cat2['photoz_best'][idx_cat2]

    total_lum=radio_luminosity(total_fluxes,z,total_SI)
    mask=total_lum > 10**24
    total_SI_loud=np.array(total_SI[mask])
    total_fluxes_loud=np.array(total_fluxes[mask])
    total_peaks_loud=np.array(total_peaks[mask])
    total_SI_quiet=np.array(total_SI[~mask])
    total_fluxes_quiet=np.array(total_fluxes[~mask])

    print('Number of quiet', len(total_fluxes_quiet))
    print('Number of loud', len(total_fluxes_loud))
    difflr,bin_centers2,median_si2,si_err2=median_SI(total_fluxes_loud,total_SI_loud,nbins=8)
    print('Spectral index loud',median_si2)
    difflr,bin_centers1,median_si1,si_err1=median_SI(total_fluxes_quiet,total_SI_quiet,nbins=5)
    print('Spectral index quiet',median_si1)

    total_fluxes=np.concatenate(fluxes)
    total_SI=np.concatenate(SI)

    # AGN_comp = (total_fluxes_loud / total_peaks_loud > 0.8) & (total_fluxes_loud / total_peaks_loud < 1.2)
    # AGN_extend = (total_fluxes_loud / total_peaks_loud > 1.5) 
    # print('Number of compact', np.sum(AGN_comp))
    # print('Number of extended', np.sum(AGN_extend))
    # difflr, bin_centers3, median_si3, si_err3 = median_SI(total_fluxes_loud[AGN_comp], total_SI_loud[AGN_comp], nbins=7)
    # difflr, bin_centers4, median_si4, si_err4 = median_SI(total_fluxes_loud[AGN_extend], total_SI_loud[AGN_extend], nbins=7)
    #total_SI, total_fluxes, bin_centers, median_si, si_err=plot(fluxes, SI)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    # Left: total_SI vs total_fluxes
    axes[0].scatter(total_fluxes, total_SI, alpha=0.5, color="#68aee0", label='All sources')  # lighter blue
    # Compute median and error for all sources
    difflr, bin_centers, median_si, si_err = median_SI(total_fluxes, total_SI, nbins=8)
    axes[0].errorbar(bin_centers, median_si, yerr=si_err, color='black', fmt='s', markersize=7, alpha=1,
                     capsize=8, linewidth=2, markeredgecolor='black', markeredgewidth=1.5, label='Median')
    axes[0].set_ylabel(r"Spectral Index $\alpha$", size=23)
    axes[0].set_xlabel("Flux density [mJy]", size=23)
    axes[0].set_ylim(-2.5, 2.5)
    axes[0].set_xlim(total_fluxes.min(), 100)
    axes[0].axhline(-1, color='red', linestyle='dotted', linewidth=3)
    axes[0].axhline(-0.3, color='red', linestyle='dotted', linewidth=3)
    axes[0].tick_params(axis='both', which='major', labelsize=17, length=8, width=1.5)
    axes[0].tick_params(axis='both', which='minor', labelsize=17, length=4, width=1)
    axes[0].set_xscale('log')
    axes[0].legend(fontsize=18)

    # Right: loud and quiet
    axes[1].scatter(total_fluxes_loud, total_SI_loud, alpha=0.8, color='orange', s=70, marker='+', label=r'$L_{1.4 \text{GHz}}> 10^{24}$ (RL)')
    axes[1].scatter(total_fluxes_quiet, total_SI_quiet, alpha=0.5, color='green', s=120, marker='*', label=r'$L_{1.4 \text{GHz}}< 10^{24}$ (RQ)')
    axes[1].errorbar(bin_centers2, median_si2, yerr=si_err2, color='orange', fmt='s', markersize=9, alpha=1,
                    capsize=8, linewidth=2, markeredgecolor='black', markeredgewidth=1.5)
    axes[1].errorbar(bin_centers1, median_si1, yerr=si_err1, color='green', fmt='s', markersize=9, alpha=1,
                    capsize=8, linewidth=2, markeredgecolor='black', markeredgewidth=1.5)
    axes[1].set_xlabel("Flux density [mJy]", size=23)
    axes[1].set_ylim(-2.5, 2.5)
    axes[1].set_xlim(total_fluxes.min(), 100)
    axes[1].axhline(-1, color='red', linestyle='dotted', linewidth=3)
    axes[1].axhline(-0.3, color='red', linestyle='dotted', linewidth=3)
    axes[1].tick_params(axis='both', which='major', labelsize=17, length=8, width=1.5)
    axes[1].tick_params(axis='both', which='minor', labelsize=17, length=4, width=1)
    axes[1].set_xscale('log')
    axes[1].legend(fontsize=18)

    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/spectral_index_flux_AGN_SFG.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()


if "__main__":
  SI=[]
  fluxes=[]
  all_ids=[]
  all_peaks=[]
  names=['A2631','Zwcl2341']
  HSC_cat='catalogs/HSC_MeerKAT_combined.fits'
  for name in names:
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/spectral/SI_cross-match.fits'
    sizes,source_ids,all_fluxes,spectral_index=find_spectral_index(real_cat)
    SI.append(spectral_index)
    fluxes.append(all_fluxes)
    all_ids.append(source_ids)
    all_peaks.append(sizes)
  plot(all_peaks,all_ids,fluxes, SI,HSC_cat)
  