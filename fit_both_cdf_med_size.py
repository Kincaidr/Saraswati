import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.interpolate import interp1d
import pickle
from astropy.table import Table

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
    theta_med_1 = median_size1(S_grid, k, m)
    theta_med_2 = median_size2(S_grid, k_fit, m_fit)

    plt.figure(figsize=(9, 7))
    plt.scatter(S_all, sizes_all, s=7, alpha=0.3, label='Resolved sources')
    plt.scatter(flux_0, sizes_0, s=15, alpha=0.3, marker='x',label='Unresolved sources',color='blue')
    plt.scatter(bin_centers, observed_medians, color='red', marker='*', s=100, label='Observed Medians')
    plt.plot(S_grid, theta_med_1, color='orange', linewidth=2, label=r"$\Theta_{\mathrm{med,1}}$")
    plt.plot(S_grid, theta_med_2, color='green', linewidth=2, label=r"$\Theta_{\mathrm{med,2}}$",)
    plt.plot(S_grid,theta_max(S_grid,sigma_low, bmaj, bmin),label="Maximum size", color="purple",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_low, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]),label="Minimum size", color="black",linestyle="--")
    plt.plot(S_grid,theta_max(S_grid,sigma_high, bmaj, bmin), color="purple",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_high, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]), color="black",linestyle="--")
    plt.xscale('log')
    plt.ylim(-5,50)
    # plt.yscale('log')
    plt.xlabel(r'$S_T$ (mJy)', size=20)
    plt.ylabel(r'$\Theta$ (arcsec)', size=20)
    plt.legend(fontsize=13)
    #plt.grid(True, which='both', ls='--')
    plt.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/size_distribution.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()

def CCDF_plot(sizes_all, bin_size=50, a=-np.log(2), b=0.62):
    theta_vals, ccdf_emp = complementary_ecdf(sizes_all, bin_size=bin_size)
    min_flux, max_flux = 0.001, 100
    S_grid = np.logspace(np.log10(min_flux), np.log10(max_flux), 1000)
    theta_lims_grid = np.logspace(np.log10(1), np.log10(100), 1000)
    
    k, m = 8, 0.3
    theta_med_1 = median_size1(S_grid, k, m)
    theta_med_2 = median_size2(S_grid, k_fit, m_fit)
    ccdf_model_1 = integral_dist(theta_med_1, theta_lims_grid, a=a, b=b)
    ccdf_model_2 = integral_dist(theta_med_2, theta_lims_grid, a=a, b=b)

    theta_lims = theta_lim(S_grid, sigma_low, bmaj, bmin, alpha_beta_array[0], alpha_beta_array[1])
    ccdf_model_1_lim = integral_dist(theta_med_1, theta_lims, a=a, b=b)
    ccdf_model_2_lim = integral_dist(theta_med_2, theta_lims, a=a, b=b)
    c1 = 1 / (1 - ccdf_model_1_lim)
    c2 = 1 / (1 - ccdf_model_2_lim)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(theta_vals, ccdf_emp, drawstyle='steps-post', linewidth=3)
    ax1.plot(theta_lims_grid, ccdf_model_1, '-', label=r"$\Theta_{\mathrm{med,1}}$", linewidth=3,color='black')
    ax1.plot(theta_lims_grid, ccdf_model_2, '-', label=r"$\Theta_{\mathrm{med,2}}$", linewidth=3,color='purple')
    ax1.axvline(k_fit, ls='--', linewidth=2, color='black')
    ax1.set_xlabel(r"$\Theta_{\mathrm{lim}} \, \mathrm{[arcsec]}$", size=18)
    ax1.set_ylabel(r"$h(>\Theta_{\mathrm{lim}})$", size=18)
    ax1.set_xlim(5, 25)
    ax1.legend(fontsize=15)
    ax1.tick_params(axis='both', which='both', direction='in', length=6, width=1)

    ax2.plot(S_grid, c1, label=r"$\Theta_{\mathrm{med,1}}$", linewidth=3,color='orange')
    ax2.plot(S_grid, c2, label=r"$\Theta_{\mathrm{med,2}}$", linewidth=3,color='green')
    ax2.set_xlabel(r"$S_T$ [mJy]", size=18)
    ax2.set_ylabel(r"$c = 1/[1 - h(>\Theta_{\mathrm{lim}})]$", size=18)
    ax2.set_xscale('log')
    ax2.legend(fontsize=15)
    ax2.tick_params(axis='both', which='both', direction='in', length=6, width=1)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolution_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.tight_layout()
    plt.show()

    interp_Func1 = interp1d(S_grid, c1, bounds_error=False, fill_value='extrapolate')
    with open('/home/kincaid/Desktop/Saraswati_codes/resolution_interp_func.pkl', 'wb') as f:
        pickle.dump(interp_Func1, f)
    return interp_Func1
    
if __name__=="__main__":
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
    sigma_low=0.016
    bmaj=8.7159409207893
    bmin=7.209874015964243
    alpha_beta_array= [1.1, 1.35]
    alpha_beta_array= [1.1, 1.7] #[1.1, 1.35] #
    bin_centers1, medians1 = equal_source_bins(flux[mask_low], sizes[mask_low], nbins=5)
    bin_centers2, medians2 = equal_source_bins(flux[~mask_low], sizes[~mask_low], nbins=6)
    bin_centers = np.concatenate([bin_centers1, bin_centers2])
    median_sizes = np.concatenate([medians1, medians2])

    # Fit
    k_fit, m_fit = fit_k_m_joint(flux, sizes,bin_centers, median_sizes, bin_size=50)
    #k_fit, m_fit=8,0.03
    # Plot
    size_dist(flux, sizes,flux_0,sizes_0,  bin_centers, median_sizes, k_fit, m_fit)
    #interp_func=CCDF_plot(sizes)
    #correction(interp_func)
    # Output
    print(f"Fitted parameters: k = {k_fit:.4f}, m = {m_fit:.4f}")
