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

    mask=total_fluxes < 100
    total_SI=np.array(total_SI[mask])
    total_fluxes=np.array(total_fluxes[mask])
    total_ids=np.array(total_ids[mask])
    common_ids,  idx_cat1, idx_cat2 = np.intersect1d(total_ids, cat2['Source_id'], return_indices=True)
    total_fluxes=total_fluxes[idx_cat1] 
    total_SI=total_SI[idx_cat1]
    z=cat2['Redshift (z)'][idx_cat2]
    total_lum=radio_luminosity(total_fluxes,z,total_SI)
    mask=total_lum > 10**23    

    total_SI_loud=np.array(total_SI[mask])
    total_fluxes_loud=np.array(total_fluxes[mask])
    total_SI_quiet=np.array(total_SI[~mask])
    total_fluxes_quiet=np.array(total_fluxes[~mask])

    
    # AGN_extend=(total_fluxes/total_peaks > 2) 
    print('Number of quiet', len(total_fluxes_quiet))
    print('Number of loud', len(total_fluxes_loud))
    difflr,bin_centers2,median_si2,si_err2=median_SI(total_fluxes_loud,total_SI_loud,nbins=6)
    print('Spectral index loud',median_si2)
    difflr,bin_centers1,median_si1,si_err1=median_SI(total_fluxes_quiet,total_SI_quiet,nbins=3)
    print('Spectral index quiet',median_si1)

    total_fluxes=np.concatenate(fluxes)
    total_SI=np.concatenate(SI)

    AGN_comp=(total_fluxes/total_peaks > 0.95) & (total_fluxes/total_peaks < 1.05) 
    print('Number of compact', np.sum(AGN_comp))
    difflr,bin_centers3,median_si3,si_err3=median_SI(total_fluxes[AGN_comp],total_SI[AGN_comp],nbins=5)
    difflr,bin_centers4,median_si4,si_err4=median_SI(total_fluxes[~AGN_comp],total_SI[~AGN_comp],nbins=5)

    plt.figure(figsize=(12, 8)) 
    plt.scatter(total_fluxes,total_SI,alpha=0.03,color='black')
    #plt.scatter(total_fluxes[AGN_extend],total_SI[AGN_extend],alpha=0.8,color='blue')
    plt.scatter(total_fluxes_loud,total_SI_loud,alpha=0.8,color='orange',s=70,marker='+',label=r'$L_{1.4 \text{GHz}}> 10^{23}$')
    plt.scatter(total_fluxes_quiet,total_SI_quiet,alpha=0.5,color='green',s=120,marker='*',label=r'$L_{1.4 \text{GHz}}< 10^{23}$')
    plt.errorbar(bin_centers2,median_si2,yerr=si_err2,color='orange',fmt='s',markersize=9,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    # plt.errorbar(bin_centers3,median_si3,yerr=si_err3,color='blue',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    # plt.errorbar(bin_centers4,median_si4,yerr=si_err4,color='purple',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    plt.errorbar(bin_centers1,median_si1,yerr=si_err1,color='green',fmt='s',markersize=9,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    plt.ylabel(r"Spectral Index $\alpha$", size=23)
    plt.xlabel("Flux density [mJy]", size=23)
    plt.ylim(-2.5,2.5)
    plt.xlim(total_fluxes.min(),100)
    plt.axhline(-1, color='red', linestyle='dotted',linewidth='3')
    plt.axhline(-0.3, color='red', linestyle='dotted',linewidth='3')
    plt.legend(fontsize=22)
    plt.tick_params(axis='both', which='major', labelsize=22, length=5, width=1)
    plt.tick_params(axis='both', which='minor', labelsize=22, length=5, width=1)
    plt.xscale('log')
    plt.savefig('plots/spectral_index_flux_AGN_SFG.png', bbox_inches='tight', pad_inches=0.1,dpi=300)   
    plt.show()

if "__main__":
  SI=[]
  fluxes=[]
  all_ids=[]
  all_peaks=[]
  names=['A2631','Zwcl2341']
  AGN_SFR_NED_cat='A2631_Zwcl2341_combined_NED_new.fits'
  for name in names:
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/spectral/SI_cross-match.fits'
    sizes,source_ids,all_fluxes,spectral_index=find_spectral_index(real_cat)
    SI.append(spectral_index)
    fluxes.append(all_fluxes)
    all_ids.append(source_ids)
    all_peaks.append(sizes)
  plot(all_peaks,all_ids,fluxes, SI,AGN_SFR_NED_cat)
  