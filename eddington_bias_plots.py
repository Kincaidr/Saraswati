import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.table import Table
from SEMPER_source_counts import SEMPER_SFG_AGN_counts
from scipy.interpolate import interp1d
from TRECS_source_counts import TRECS_counts
from matplotlib.lines import Line2D

def poly_fit(x,y):
    coeffs = np.polyfit(x, y, 6)
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

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def gamma_from_SCs_Bondi(S, a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a = np.array([a0, a1, a2, a3, a4, a5, a6])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def gamma_from_SCs(S, coef):
    # a = np.array([coef[6], coef[6], coef[5], coef[3], coef[4], coef[5], coef[6],coef[7]])
    a = coef[::-1]
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def compute_eddington_bias(name,line,gamma_values):
    SNR_values=[5,8,20]
    flux_array=np.linspace(0.01,100,10000) # Convert to mJy    
    f1=interp1d(S_values_gamma, gamma_values, bounds_error=False, fill_value="extrapolate")
    ax1,ax2 = axes
    ax1.plot(S_values_gamma, gamma_values,label=name, color='black',linestyle=line)
    ax1.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    ax1.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('linear')
    ax1.set_ylabel(r'$ \gamma$',size=20 )
    ax1.set_ylim(1.5, 2.6)
    ax1.axvline(11e-3, color='purple', linestyle='-')
    ax1.axvline(16e-3, color='orange', linestyle='-')
    #ax1.set_xlabel('Flux Density (mJy)',size=17)
    ax1.legend(fontsize=16, loc='upper right')
    ax2 = axes[1]
    for SNR in SNR_values:
        S_true=(flux_array/2)*(1+np.sqrt(1-(4* f1(flux_array)/(SNR**2))))
        y = flux_array / S_true
        ax2.plot(flux_array,y,label=f'SNR={SNR}',alpha=1,linestyle=line,color='black')

        if name == 'Bondi +08':
            idx = np.abs(flux_array - 1).argmin()
            ax2.text(flux_array[idx], y[idx] + 0.01, f'SNR={SNR}', color='purple',
                    fontsize=15, ha='left', va='center')
    ax2.set_xscale('log')
    ax2.set_yscale('linear')
    ax2.set_ylabel(r'$S_{obs}/S_{true}$',size=20)
    ax2.set_xlabel('Flux Density (mJy)',size=20 )
    ax2.axvline(11e-3, color='purple', linestyle='-')
    ax2.axvline(16e-3, color='orange', linestyle='-')     
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)   
    #ax2.set_xlim(0.1, 100)
    ax2.set_ylim(0.95, 1.16)
    

# def plot_gamma(coef):
#     S_values_gamma = np.linspace(0.1, 100, 1000)  # Flux values from 0.001 to 1000 mJy
#     gamma_values_Semper = gamma_from_SCs_SEMPER(S_values_gamma,coef)
#     gamma_values_Bondi = gamma_from_SCs_Bondi(S_values_gamma)
#     f1=interp1d(S_values, gamma_values_Bondi, bounds_error=False, fill_value="extrapolate")
#     f2=interp1d(S_values, gamma_values_Semper, bounds_error=False, fill_value="extrapolate")
#     plt.plot(S_values, gamma_values_Bondi, label='Bondi fit', color='black',linestyle='--')
#     plt.plot(S_values, gamma_values_Semper,label='SEMPER fit', color='black',linestyle='-')
#     plt.xscale('log')
#     plt.yscale('linear')
#     plt.ylabel(r'$ \gamma$')
#     plt.xlabel('Flux Density (mJy)',size=15)
#     plt.legend()
#     plt.grid(True)
#     plt.show()
#     return(f1,f2)

if __name__ == "__main__":
    outpath='/home/kincaid/Desktop/Saraswati_codes/'
    names=['Bondi +08','Semper +25','TRECS +23']
    linestyles = ['-', '--', '-.']
    
    fig, axes = plt.subplots(2,1, figsize=(10, 8))  # Create a side-by-side layout
    for (i,name),line in zip(enumerate(names),linestyles):
        if name == 'Bondi +08':
            S_values_gamma = np.linspace(0.01, 100, 10000)
            gamma_values = gamma_from_SCs_Bondi(S_values_gamma)
        elif name == 'Semper +25':
            x, y = SEMPER_SFG_AGN_counts()
            mask = (x < 300) & (y >0) 
            x = np.log10(x[mask])
            y = np.log10(y[mask])
            coef=poly_fit(x,y)
            gamma_values = gamma_from_SCs(S_values_gamma,coef)
        elif name == 'TRECS +23':
            x, y = TRECS_counts()
            mask = (x < 300) & (y >0)
            x = np.log10(x[mask])
            y = np.log10(y[mask])
            coef=poly_fit(x,y)
            gamma_values = gamma_from_SCs(S_values_gamma,coef)
    #interp_gamma_bondi,interp_gamma_semper=plot_gamma(coef)
        compute_eddington_bias(name,line,gamma_values)
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/eddington_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
    