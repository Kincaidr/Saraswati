from astropy.table import Table
from scipy.stats import linregress
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def find_spectral_index(real_cat):
    cat = Table.read(real_cat)
    freqs = [998, 1283, 1569]  # Using the three frequencies

    spectral_index = []
    for i in range(len(cat['Total_flux'])):
        all_fluxes = [cat['Total_flux_1'][i], cat['Total_flux_2'][i], cat['Total_flux'][i]]
        all_freqs = [freqs[0], freqs[1], freqs[2]]
        Result = linregress(np.log10(all_freqs), np.log10(all_fluxes))
        spectral_index.append(Result[0])
    return(spectral_index)

def plot(SI):
    plt.figure(figsize=(12, 8))
    plt.figure(figsize=(10, 8))
    plt.hist(SI[1], bins=30, color='red', alpha=1, density=False)
    plt.hist(SI[0], bins=30, color='blue', alpha=0.5, density=False)
    plt.xlabel(r"Spectral Index $\alpha$ ", size=23)
    plt.ylabel("Number of Sources", size=23)
    plt.axvline(np.median(SI[0]), color='blue', linestyle='dashed', label=rf'A2631, $\alpha \sim $ {np.median(SI[0]):.2f}',linewidth=2)
    plt.axvline(np.median(SI[1]), color='red', linestyle='dashed', label=rf'Zwcl2341, $\alpha \sim $ {np.median(SI[1]):.2f}',linewidth=2)
    plt.legend(fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18, length=5, width=1)
    plt.tick_params(axis='both', which='minor', labelsize=18, length=5, width=1)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/Spectral_Index_Distribution.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()


if "__main__":
  SI=[]
  names=['A2631','Zwcl2341']
  for i,name in enumerate(names):
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/spectral/SI_cross-match_new.fits'
    spectral_index=find_spectral_index(real_cat)
    SI.append(spectral_index)
    print('Number of sources for '+ name,len(SI[i]))
  plot( SI)
  