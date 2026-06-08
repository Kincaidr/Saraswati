import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from astropy.table import Table

# --- Step 1: Your existing functions ---
def median_size2(S, k, m):
    S = np.asarray(S)
    return np.where(S < 1, k, k * S**m)

def integral_dist(theta_med, theta_lim, a=-np.log(2), b=0.62):
    return np.exp(a * (theta_lim / theta_med)**b)

def complementary_ecdf(data, bin_size):
    counts, bin_edges = np.histogram(data, bins=bin_size, density=False)
    cumsum = np.cumsum(counts[::-1])[::-1]
    ccdf = cumsum / cumsum[0]
    return bin_edges[:-1], ccdf  # bin centers (left edges), ccdf

# --- Step 2: Fitting function ---
def fit_k_m_to_data(observed_sizes, S, bin_size=30):
    theta_lims, ccdf_empirical = complementary_ecdf(observed_sizes, bin_size)

    def loss(params):
        k, m = params
        theta_med = median_size2(S, k, m)
        if np.any(theta_med <= 0):
            return np.inf  # prevent invalid values

        theta_med_scalar = np.mean(theta_med)
        model_ccdf = integral_dist(theta_med_scalar, theta_lims)
        return np.sum((ccdf_empirical - model_ccdf)**2)

    result = minimize(loss, x0=[1.0, 5.0], bounds=[(1e-5, None), (0.01, 5)])
    k_fit, m_fit = result.x
    print(k_fit,m_fit)
    # Plotting
    theta_med_fit = np.mean(median_size2(S, k_fit, m_fit))
    model_ccdf = integral_dist(theta_med_fit, theta_lims)
    breakpoint()
    plt.figure()
    plt.plot(theta_lims, ccdf_empirical, drawstyle='steps-post', label='Empirical CCDF')
    plt.plot(theta_lims, model_ccdf, '-', label=f'Model Fit: k={k_fit:.3f}, m={m_fit:.3f}')
    plt.xlabel('Observed Source Size (θ)')
    plt.ylabel('CCDF')
    plt.xlim(0,50)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return k_fit, m_fit

def info(cat, min_flux, max_flux):
    mask=(cat['Maj']==0)
    flux=cat['Total_flux'][~mask]*1e3
    min=cat['Min'][~mask]*3600
    maj=cat['Maj'][~mask]*3600
    size=np.sqrt(maj*min)
    mask= (flux>min_flux) & (flux < max_flux)
    fluxes=flux[mask]
    true_sizes=size[mask]
    return(true_sizes,fluxes)

if __name__=="__main__":
    name='A2631'
    min_flux=1
    max_flux=1000
    S=np.logspace(-3,3,10000)
    cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_eddington_corr_srl.fits')
    observed_sizes=info(cat, min_flux, max_flux)
    fit_k_m_to_data(observed_sizes, S, bin_size=100)