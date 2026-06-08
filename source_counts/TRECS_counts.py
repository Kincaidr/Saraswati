from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pickle

def source_counts():
    cat='catalogue_continuum_wrapped.fits'
    table = Table.read(cat)
    flux_1400=table['I1400']
    area_sr = 25 * (np.pi / 180) ** 2
    flux=np.log10(flux_1400/1000)
    f,bin_edges=np.histogram(flux, bins=40,range=(-7,4))
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    s=10**bin_centers
    bin_size = bin_edges[1] - bin_edges[0]
    factor = 0.434294482  # 1/ln(10)
    dnds = factor / s * f / bin_size  # dN/dlog10S
    f14_s = (dnds * s**2.5) / area_sr  # S^2.5 dN/dS
    f_err = np.sqrt(f)
    dnds_up = factor / s * (f + f_err) / bin_size
    dnds_down = factor / s * (f - f_err) / bin_size

    f14_s_up = (dnds_up * s**2.5) / area_sr
    f14_s_down = (dnds_down * s**2.5) / area_sr
    TRECS_func = interp1d(bin_centers, f14_s, bounds_error=False, fill_value=0)
    with open('TRECS_func.pkl', 'wb') as f:
        pickle.dump(TRECS_func, f)
    return (10**bin_centers)*1e3, f14_s
    # plt.errorbar((10**bin_centers)*1e3, f14_s, 
    #             yerr=[f14_s - f14_s_down, f14_s_up - f14_s], 
    #             fmt='o')
    # plt.scatter((10**bin_centers)*1e3, f14_s, color='blue')
    # plt.xlabel('Flux [mJy]')
    # plt.ylabel(r'$S^{2.5} \frac{dN}{dS}$ [sr$^{-1}$ Jy$^{1.5}$]')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.show()

if  __name__ == "__main__":
    source_counts()