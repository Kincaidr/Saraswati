import bdsf
import json
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)


def ratio_source_count(flux, nbins, equal_source_size=False, equal_bin_size=False, Range_x= None):

        if  equal_source_size:
                pers = np.linspace(0, 100, nbins+1)
                Range_x = np.percentile(flux, pers)
        elif equal_bin_size:
                Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1)
        else:
              Range_x=Range_x
        centres = (Range_x[:-1] + Range_x[1:]) / 2.0
        hist, _ = np.histogram(flux, bins=Range_x)
        source_tot=hist
        return(Range_x, source_tot, centres)
    

def ratio(merged_rec,merged_inv,plots, nbins, noise, output):
   
    rec_cat = Table.read(merged_rec)
    inv_cat = Table.read(merged_inv)
    flux_rec=rec_cat['Total_flux']
    flux_inv=inv_cat['Total_flux'][inv_cat['Total_flux']> 3*noise*10**-3]
    Range_x,counts_inv, bin_centres=ratio_source_count(flux_inv, nbins=nbins, equal_source_size=True, equal_bin_size=False, Range_x= None)
    _,counts_rec, _=ratio_source_count(flux_rec, nbins=nbins, equal_source_size=False, equal_bin_size=False, Range_x= Range_x)
    
    ratio=counts_inv/counts_rec
    mask=ratio< 1
    ratio=ratio[mask]
 
    #error=ratio*np.sqrt((np.sqrt(counts_inv)/counts_inv)**2+(np.sqrt(counts_rec)/counts_rec)**2)
    print(ratio)
    #plt.errorbar(bin_centres*10**3,ratio,yerr=error) 
    breakpoint()
    
    bins=np.array(len(counts_inv))
 
    plt.hist(np.arange(len(counts_inv)), weights=counts_inv, bins=bins, alpha=0.5, color='blue', label='counts_inv', edgecolor='black', density=False)
    plt.hist(np.arange(len(counts_rec)), weights=counts_rec, bins=bins, alpha=0.5, color='gray', label='counts_rec', edgecolor='black', density=False)
    plt.xlabel('Log Flux S [mJy]',size=12)
    plt.ylabel('Unormalized Counts',size=12)
    plt.legend()
    plt.show()
    plt.savefig(output+"false_detections_unormalized.png")
    print(plots+'false_detections_unormalized.png saved')
    plt.close()
  

    plt.hist(x=flux_rec*10**3,bins=Range_x*10**3, alpha=0.5, color='gray', label='counts_rec', edgecolor='black', density=False)
    plt.hist(x=flux_inv*10**3, bins=Range_x*10**3, alpha=0.8, color='blue', label='counts_inv', edgecolor='black', density=False)
    #plt.xticks(ticks=np.arange(len(Range_x)), labels=Range_x, rotation=0)
    plt.xlabel('Log Flux S [mJy]',size=12)
    plt.ylabel('normalized Counts',size=12)
    plt.legend()
    plt.show()
    plt.savefig(output+"false_detections_normalized.png")
    print(output+'false_detections_unormalized.png saved')
    plt.close()


if __name__ == "__main__":

    config_path = 'A2631.json'
    config = read_config(config_path)
    path=config['path']
    noise=config['sigma']
    catalogs=config['catalogs']
    output=path+catalogs
    merged_rec='/home/kincaid/Desktop/Saraswati_codes/A2631/catalogs/merged_rec.fits'
    merged_inv='/home/kincaid/Desktop/Saraswati_codes/A2631/catalogs/merged_inverted.fits'
    plots=config['plots']
    nbins=20
    #new_fits_image=invert_image(fits_image)
    #cat=catalog_generation(new_fits_image)
    ratio(merged_rec,merged_inv,plots, nbins,noise, output)
    #analysis(cat, name,noise)