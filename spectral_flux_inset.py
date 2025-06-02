from astropy.table import Table
from scipy.stats import linregress
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
  
    mask=total_fluxes < 100
    breakpoint()
    nbins=10
    difflr,bin_centers,median_si,si_err=median_SI(total_fluxes[mask],total_SI[mask],nbins)
    plt.figure(figsize=(12, 8)) 
    plt.scatter(total_fluxes[mask],total_SI[mask],alpha=0.5)
    plt.errorbar(bin_centers,median_si,yerr=si_err,color='black',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5 )
    plt.ylabel(r"Spectral Index $\alpha$ ", size=23)
    plt.xlabel("Flux density [mJy]", size=23)
    plt.ylim(-3,3)
    plt.axhline(-1, color='red', linestyle='dotted',linewidth='3')
    plt.axhline(-0.3, color='red', linestyle='dotted',linewidth='3')
    plt.legend(fontsize=17)
    plt.tick_params(axis='both', which='major', labelsize=22, length=5, width=1)
    plt.tick_params(axis='both', which='minor', labelsize=22, length=5, width=1)
    plt.xscale('log')

    ax_inset = inset_axes(plt.gca(), width="45%", height="40%", loc='upper right', borderpad=2)
    total_lum=radio_luminosity(total_fluxes,0.27,total_SI)
    mask=total_lum > 10**23
    total_SI_loud=np.array(total_SI[mask])
    total_fluxes_loud=np.array(total_fluxes[mask])
    total_SI_quiet=np.array(total_SI[~mask])
    total_fluxes_quiet=np.array(total_fluxes[~mask])
    difflr,bin_centers2,median_si2,si_err2=median_SI(total_fluxes_loud,total_SI_loud,nbins=6)
    difflr,bin_centers1,median_si1,si_err1=median_SI(total_fluxes_quiet,total_SI_quiet,nbins=3)

    ax_inset.scatter(total_fluxes,total_SI,alpha=0.03,color='black')

    #plt.scatter(total_fluxes[AGN_comp],total_SI[AGN_comp],alpha=0.8,color='black')
    #plt.scatter(total_fluxes[AGN_extend],total_SI[AGN_extend],alpha=0.8,color='blue')
    ax_inset.scatter(total_fluxes_quiet,total_SI_quiet,alpha=0.5,color='green',s=20,marker='*',label='Radio quiet')
    ax_inset.scatter(total_fluxes_loud,total_SI_loud,alpha=0.8,color='orange',s=20,marker='+',label='Radio loud')
    ax_inset.errorbar(bin_centers2,median_si2,yerr=si_err2,color='orange',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    #plt.errorbar(bin_centers3,median_si3,yerr=si_err3,color='blue',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    ax_inset.errorbar(bin_centers1,median_si1,yerr=si_err1,color='green',fmt='s',markersize=7,alpha=1 ,capsize=8, linewidth=2, markeredgecolor='black',markeredgewidth=1.5)
    #ax_inset.set_ylabel(r"Spectral Index $\alpha$", size=23)
    #ax_inset.set_xlabel("Flux density [mJy]", size=23)
    ax_inset.set_ylim(-3,3)
    ax_inset.axhline(-1, color='red', linestyle='dotted',linewidth='3')
    ax_inset.axhline(-0.3, color='red', linestyle='dotted',linewidth='3')
    ax_inset.legend(fontsize=13)
    #ax_inset.tick_params(axis='both', which='major', labelsize=12, length=5, width=1)
    #ax_inset.tick_params(axis='both', which='minor', labelsize=12, length=5, width=1)
    ax_inset.set_xscale('log')
    plt.savefig('plots/spectral_index_flux_inset.png', bbox_inches='tight', pad_inches=0.1)
    plt.show()

if "__main__":
  SI=[]
  fluxes=[]
  names=['A2631','Zwcl2341']
  for name in names:
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/spectral/SI_cross-match.fits'

    all_fluxes,spectral_index=find_spectral_index(real_cat)
    SI.append(spectral_index)
    fluxes.append(all_fluxes)
  plot(fluxes, SI)