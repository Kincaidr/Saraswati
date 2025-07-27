import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.interpolate import interp1d
import pickle

def func(a,b,S,N):
    f=a+b/(S/N)
    return(f)

def median_size1(S, k, m):
    t = k * S**m
    return t

def theta_max(S,sigma, bmaj, bmin):
    theta_max=np.sqrt(bmaj*bmin)*np.sqrt((S/(threshold*sigma)))
    return(theta_max)

def theta_min(S,N, bmaj, bmin,a,b):
    mean_beam=np.sqrt(bmaj*bmin)
    theta_min=mean_beam*np.sqrt(func(a,b,S,N)-1)
    return(theta_min)

def theta_lim(S, sigma, bmaj, bmin, a,b):
        x=theta_min(S, sigma, bmaj, bmin, a,b)
        y=theta_max(S, sigma, bmaj, bmin)
        theta_lim=np.maximum(x,y)
        return(theta_lim)

def median_size2(S, k, m):
    """Median size model as a function of flux S."""
    S = np.asarray(S)
    return np.where(S < 1, k, k * S**m)

def integral_dist(theta_med, theta, a=-np.log(2), b=0.62):
    """Integral distribution model (CCDF)."""
    return np.exp(a * (theta / theta_med)**b)

def complementary_ecdf(data, bin_size):
    """Empirical complementary CDF (CCDF) from histogram."""
    counts, bin_edges = np.histogram(data, bins=bin_size, density=False)
    cumsum = np.cumsum(counts[::-1])[::-1]
    ccdf = cumsum / cumsum[0]
    return bin_edges[:-1], ccdf  # bin edges (left), CCDF

def equal_source_bins(x, y, nbins=6):
    """Bin data by equal number of sources per bin."""
    x = np.asarray(x)
    y = np.asarray(y)
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    y = y[sorted_idx]
    bins = np.array_split(np.arange(len(x)), nbins)

    bin_centers = []
    median_values = []
    for b in bins:
        bin_centers.append(np.median(x[b]))
        median_values.append(np.median(y[b]))
    return np.array(bin_centers), np.array(median_values)


def total_loss(params, S_all, sizes_all, bin_centers, observed_medians, bin_size=30, w_ccdf=1.0, w_median=1.0):
    """Combined loss: CCDF fit + Median vs flux fit."""
    k, m = params
    if k <= 0 or m < 0:
        return np.inf  # Reject unphysical values

    theta_meds_all = median_size2(S_all, k, m)
    theta_med_scalar = np.mean(theta_meds_all)  # global median
    theta_vals, ccdf_emp = complementary_ecdf(sizes_all, bin_size=bin_size)
    model_ccdf = integral_dist(theta_med_scalar, theta_vals)
    loss_ccdf = np.sum((ccdf_emp - model_ccdf)**2)
    model_medians = median_size2(bin_centers, k, m)
    loss_median = np.sum((observed_medians - model_medians)**2)
    return w_ccdf * loss_ccdf + w_median * loss_median

def fit_k_m_joint(S_all, sizes_all, bin_centers, observed_medians, bin_size=30):
    """Fit k, m jointly to both CCDF and median–flux trend."""
    result = minimize(
        total_loss,
        x0=[8, 0.05],
        args=(S_all, sizes_all, bin_centers, observed_medians, bin_size),
        bounds=[(1e-4, 10), (0.01, 5)],
        method='L-BFGS-B'
    )
    return result.x  # k, m


def size_dist(S_all, sizes_all,flux_0,sizes_0, bin_centers, observed_medians, k_fit, m_fit, bin_size=30):
    S_grid = np.logspace(np.log10(np.min(S_all)), np.log10(np.max(S_all)), 1000)
    k,m=8,0.3
    theta_med_1 = median_size2(S_grid, k, m)
    theta_med_2 = median_size2(S_grid, k_fit, m_fit)

    plt.figure(figsize=(8, 5))
    plt.scatter(S_all, sizes_all, s=5, alpha=0.3, label='Resolved sources')
    plt.scatter(flux_0, sizes_0, s=15, alpha=0.3, label='Unresolved sources',color='blue')
    plt.scatter(bin_centers, observed_medians, color='red', marker='*', s=100, label='Observed Medians')
    plt.plot(S_grid, theta_med_1, color='black', linewidth=2, label=r"$\Theta_{\mathrm{med,1}}$")
    plt.plot(S_grid, theta_med_2, color='purple', linewidth=2, label=r"$\Theta_{\mathrm{med,2}}$",)
    
    plt.plot(S_grid,theta_max(S_grid,sigma_low, bmaj, bmin),label="Maximum size", color="green",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_low, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]),label="Minimum size", color="orange",linestyle="--")
    plt.plot(S_grid,theta_max(S_grid,sigma_high, bmaj, bmin), color="green",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_high, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]), color="orange",linestyle="--")
    plt.xscale('log')
    plt.ylim(0,30)
    # plt.yscale('log')
    plt.xlabel('Flux Density (mJy)')
    plt.ylabel('Source Size (arcsec)')
    plt.legend()
    #plt.grid(True, which='both', ls='--')
    plt.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    plt.tight_layout()
    plt.show()

def CCDF_plot(S_all,sizes_all,bin_size=30):
    theta_vals, ccdf_emp = complementary_ecdf(sizes_all, bin_size=bin_size)
    theta_med_2 = np.mean(median_size2(S_all, k_fit, m_fit))
    theta_med_1 = np.mean(median_size1(S_all, k_fit, m_fit))
    ccdf_model_1 = integral_dist(theta_med_1, theta_vals)
    ccdf_model_2 = integral_dist(theta_med_2, theta_vals)
    plt.figure(figsize=(8, 5))
    plt.plot(theta_vals, ccdf_emp, drawstyle='steps-post', label='Empirical CCDF')
    plt.plot(theta_vals, ccdf_model_2, '-', label=r"$\Theta_{\mathrm{med,2}}$", linewidth=2)
    plt.plot(theta_vals, ccdf_model_1, '-', label=r"$\Theta_{\mathrm{med,1}}$", linewidth=2)
    # plt.xscale('log')
    # plt.yscale('log')
    #   plt.xlim(0,30)
    plt.xlabel('Source Size (arcsec)')
    plt.ylabel('CCDF')
    plt.legend()
    plt.grid(True, which='both', ls='--')
    plt.tight_layout()
    plt.show()
    
def resol_bias():
    S_vals = np.logspace(-2, 2, 1000) 
    theta_lims = theta_lim(S_vals, sigma_low, bmaj, bmin, alpha_beta_array[0],alpha_beta_array[1])
    theta_med_1 = median_size2(S_vals, k_fit, m_fit)
    theta_med_2 = median_size1(S_vals, k_fit, m_fit)
    h1=integral_dist(theta_med_1,theta_lims)
    h2=integral_dist(theta_med_2,theta_lims)
    c1=1/(1-h1)
    c2=1/(1-h2)
    interp_Func1 = interp1d(S_vals, c2, bounds_error=False, fill_value=0)
    with open('/home/kincaid/Desktop/Saraswati_codes/resolution_interp_func.pkl', 'wb') as f:
        pickle.dump(interp_Func1, f)
    plt.plot(S_vals ,c1,label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    plt.plot(S_vals ,c2,label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    plt.ylabel(r"c=$1/[1-h(>\Theta_{\mathrm{lim}})$]", size=24)
    plt.xlabel(r"$S_T$ [mJy]", size=24)
    plt.xscale('log')
    plt.legend(fontsize=22)
    plt.show()
    return interp_Func1

def correction(interp_func):
    mask=cat['Maj'] !=0
    flux=cat['Total_flux'][mask]*1e3
    corr=interp_func(np.array(flux))
    S_corr=flux*1/corr
    cat['Total_flux'][mask] = S_corr*1e-3  # Convert back to Jy
    outname=f'{path}{name}/catalogs/{name}_resolution_corr_srl.fits'
    cat.write(outname, overwrite=True)
    print(f"Resolution bias corrected catalog written to {outname}")

from astropy.table import Table

if __name__=="__main__":
    name='Zwcl2341'
    path='/home/kincaid/Desktop/Saraswati_codes/'
    cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/MeerKAT_eddington_combined.fits')
    mask = cat['Maj'] == 0
    flux = cat['Total_flux'][~mask] * 1e3  # mJy
    minor = cat['Min'][~mask] * 3600  # arcsec
    major = cat['Maj'][~mask] * 3600  # arcsec
    sizes = np.sqrt(major * minor)    # geometric mean size
    flux_0 = cat['Total_flux'][mask] * 1e3  # mJy
    minor_0 = cat['Min'][mask] * 3600  # arcsec
    major_0 = cat['Maj'][mask] * 3600  # arcsec
    sizes_0 = np.sqrt(major_0 * minor_0)

    # Bin flux for observed medians
    threshold=5
    mask_low = flux < 1
    sigma_high=0.04
    sigma_low=0.015
    bmaj=8.7159409207893
    bmin=7.209874015964243
    alpha_beta_array= [1.1, 1.7] #[1.1, 1.35] #
    bin_centers1, medians1 = equal_source_bins(flux[mask_low], sizes[mask_low], nbins=5)
    bin_centers2, medians2 = equal_source_bins(flux[~mask_low], sizes[~mask_low], nbins=6)
    bin_centers = np.concatenate([bin_centers1, bin_centers2])
    median_sizes = np.concatenate([medians1, medians2])

    # Fit
    k_fit, m_fit = fit_k_m_joint(flux, sizes,bin_centers, median_sizes, bin_size=50)
    k_fit, m_fit=8,0.03
    # Plot
    size_dist(flux, sizes,flux_0,sizes_0,  bin_centers, median_sizes, k_fit, m_fit)
    CCDF_plot(flux,sizes,bin_size=30)
    interp_func=resol_bias()
    correction(interp_func)
    # Output
    print(f"Fitted parameters: k = {k_fit:.4f}, m = {m_fit:.4f}")
