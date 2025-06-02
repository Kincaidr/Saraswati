import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt


def plot(cat1,cat2):
    flux1=Table.read(cat1)
    flux2=Table.read(cat2)
    flux1=flux1['Total_flux_inj']
    flux2=flux2['Total_flux']
    breakpoint()
    plt.scatter(flux1,flux2)
    plt.plot(flux1,flux1)
    plt.xlabel("inj flux")
    plt.ylabel("rec flux")
    plt.show()
    plt.savefig('Rec_vs_inj_flux_pybdsf.png')

if __name__=="__main__":
 
    inj_cat='injected_cat_n0.fits'
    rec_cat='simulated_image_test_n0.fits_srl.fits'

    plot(inj_cat,rec_cat)
    










   
 
 
 
 


