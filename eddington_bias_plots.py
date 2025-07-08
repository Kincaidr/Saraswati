import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.table import Table
from SEMPER_source_counts import SEMPER_SFG_AGN_counts
from scipy.interpolate import interp1d
from TRECS_source_counts import TRECS_counts


def poly_fit():
    Semper_M,Semper_counts=  SEMPER_SFG_AGN_counts()
    x=Semper_M
    y=np.log10(Semper_counts)
    breakpoint()
    coeffs = np.polyfit(x, y, 7)
    p = np.poly1d(coeffs)
    print("Coefficients of the 7th order polynomial fit:", coeffs)
    x_fit = np.linspace(min(x), max(x), 1000)
    y_fit = p(x_fit)
    # plt.scatter(x, y, label='Data')
    # plt.plot(x_fit, y_fit, 'r-', label='7th Order Fit')
    # plt.legend()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('7th Order Polynomial Fit')
    # plt.grid(True)
    # plt.show()
    return coeffs

def gamma_from_SCs(S, a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a = np.array([a0, a1, a2, a3, a4, a5, a6])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def gamma_from_SCs_SEMPER(S, coef):
    a = np.array([coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6],coef[7]])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def compute_eddington_bias(interp_gamma):
    SNR_values=[5,10,20]
    flux_array=np.linspace(0.01,100,1000) # Convert to mJy
    for SNR in SNR_values:
        S_true=(flux_array/2)*(1+np.sqrt(1-(4*interp_gamma(flux_array)/(SNR**2))))
        plt.plot(flux_array,flux_array/S_true,label=f'SNR={SNR}',alpha=0.5)
    plt.xscale('log')
    plt.yscale('linear')
    plt.ylabel(r'$S_{obs}/S_{true}$',size=17)
    plt.xlabel('Flux Density (mJy)',size=17)
    plt.xlim(0.01, 100)
    plt.ylim(0.6, 1.2)
    plt.legend()
    plt.show()
    

def plot_gamma(coef):
    S_values = np.linspace(0.1, 100, 1000)  # Flux values from 0.001 to 1000 mJy
    gamma_values = gamma_from_SCs_SEMPER(S_values,coef)
    #gamma_values = gamma_from_SCs(S_values)
    f=interp1d(S_values, gamma_values, bounds_error=False, fill_value="extrapolate")
    plt.plot(S_values, gamma_values)
    plt.xscale('log')
    plt.yscale('linear')
    plt.ylabel(r'$ \gamma$')
    plt.xlabel('Gamma')
    plt.legend()
    plt.grid(True)
    plt.show()
    return(f)

if __name__ == "__main__":
    path='/home/kincaid/Desktop/Saraswati_codes/'
    coef=poly_fit()
    interp_gamma=plot_gamma(coef)
    compute_eddington_bias(interp_gamma)
