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
    return L

def find_spectral_index(real_cat):
    cat = Table.read(real_cat)
    freqs = [998, 1283, 1569]  # Using the three frequencies
    spectral_index = []
    fluxes=[]
    for i in range(len(cat['Total_flux'])):
        all_fluxes = [cat['Total_flux_1'][i], cat['Total_flux_2'][i], cat['Total_flux'][i]]
        all_freqs = [freqs[0], freqs[1], freqs[2]]
        Result = linregress(np.log10(all_freqs), np.log10(all_fluxes))
        spectral_index.append(Result[0])
        fluxes.append(cat['Total_flux_2'][i]*1e3)
    return(fluxes,spectral_index)

def median_SI(fluxes, SI,nbins):
    Range_x = 10**np.linspace(start=np.log10(fluxes.min()), stop=np.log10(fluxes.max()), num=nbins+1)
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

def plot(fluxes,SI):
    total_SI=np.concatenate(SI)
    total_fluxes=np.concatenate(fluxes)
    mask=total_fluxes < 80
    total_SI=np.array(total_SI[mask])
    total_fluxes=np.array(total_fluxes[mask])
    nbins=10
    difflr,bin_centers,median_si,si_err=median_SI(total_fluxes,total_SI,nbins)
    plt.figure(figsize=(8, 6)) 
    plt.scatter(total_fluxes,total_SI,alpha=0.5)
    plt.errorbar(bin_centers,median_si,yerr=si_err,color='black',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5 )
    plt.ylabel(r"Spectral Index $\alpha$ ", size=23)
    plt.xlabel("Flux density [mJy]", size=23)
    plt.ylim(-2.5,2.5)
    plt.axhline(-1, color='red', linestyle='dotted',linewidth='3')
    plt.axhline(-0.3, color='red', linestyle='dotted',linewidth='3')
    plt.legend(fontsize=17)
    plt.tick_params(axis='both', which='major', labelsize=22, length=5, width=1)
    plt.tick_params(axis='both', which='minor', labelsize=22, length=5, width=1)
    plt.xscale('log')
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/spectral_index_flux.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
    return total_SI, total_fluxes, bin_centers, median_si, si_err

if "__main__":
  SI=[]
  fluxes=[]
  names=['A2631','Zwcl2341']
  for name in names:
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/spectral/SI_cross-match_new.fits'
    all_fluxes,spectral_index=find_spectral_index(real_cat)
    SI.append(spectral_index)
    fluxes.append(all_fluxes)
  total_SI, total_fluxes, bin_centers, median_si, si_err=plot(fluxes, SI)