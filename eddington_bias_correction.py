import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.table import Table
from SEMPER_source_counts import SEMPER_SFG_AGN_counts
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D

def poly_fit():
    Semper_M,Semper_counts=  SEMPER_SFG_AGN_counts()
    x=np.log10(Semper_M)
    y=np.log10(Semper_counts)
    coeffs = np.polyfit(x, y, 7)
    p = np.poly1d(coeffs)
    print("Coefficients of the 7th order polynomial fit:", coeffs)
    x_fit = np.linspace(min(x), max(x), 1000)
    y_fit = p(x_fit)
    plt.scatter(x, y, label='Data')
    plt.plot(x_fit, y_fit, 'r-', label='7th Order Fit')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('7th Order Polynomial Fit')
    plt.grid(True)
    plt.show()
    return coeffs

def gamma_from_SCs_Bondi(S, a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a = np.array([a0, a1, a2, a3, a4, a5, a6])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def gamma_from_SCs_SEMPER(S, coef):
    # a = np.array([coef[6], coef[6], coef[5], coef[3], coef[4], coef[5], coef[6],coef[7]])
    a = coef[::-1]
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def compute_eddington_bias(catalog):
    cat=Table.read(catalog)
    mask=cat['Total_flux'] ==cat['Peak_flux']
    print(f"Fraction of unresolved sources that willl be corrected: {np.sum(mask)/len(cat)*100:.2f}%")
    flux = cat['Total_flux'][mask] 
    rms=cat['Isl_rms'][mask]  
    SNR=flux/rms
    S_true=(flux/2)*(1+np.sqrt(1-(4*gamma_from_SCs_SEMPER(flux, coef)/(SNR**2))))
    #S_true=(flux/2)*(1+np.sqrt(1-(4*gamma_from_SCs_Bondi(flux)/(SNR**2))))
    plt.hist(np.log10(cat['Total_flux']), bins=100, alpha=0.7, label='Total Flux before',color='blue')
    cat['Total_flux'][mask] = S_true
    outname=f'{path}{name}/catalogs/{name}_eddington_corr_srl.fits'
    cat.write(outname, overwrite=True)
    print(f"Eddington bias corrected catalog written to {outname}")
    cat=Table.read(outname)
    plt.hist(np.log10(cat['Total_flux']), bins=100, alpha=0.7, label='Total Flux after',color='red')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    path='/home/kincaid/Desktop/Saraswati_codes/'
    coef=poly_fit()
    names=['A2631','Zwcl2341']

    for name in names:
        catalog= f'{path}/{name}/catalogs/{name}_flux_corr_srl.fits'
        compute_eddington_bias(catalog)
