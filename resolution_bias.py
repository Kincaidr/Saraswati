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

def h_bondi(phi):                  
    return (1./(1.6)**phi)*(phi <= 4.) +  (phi**(-1.3)-0.01)*(phi>4.)

def median_size2(S, k, m):
    S = np.asarray(S)
    return np.where(S < 1, k, k * S**m)


# def integral_dist(theta_med, theta, a=-np.log(2), b=0.62):
#     return np.exp(a * (theta / theta_med)**b)

def integral_dist(theta_med, theta, a, b):
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


def loss(params, S, theta_obs):
    k, m = params
    theta_model = median_size2(S, k, m)
    return np.sum((theta_model - theta_obs)**2)


def ccdf_loss(params, theta_vals, ccdf_emp, theta_med):
    a, b = params
    model = integral_dist(theta_med, theta_vals, a, b)
    mask = (ccdf_emp > 0) & (model > 0)
    return np.sum((np.log10(ccdf_emp[mask]) - np.log10(model[mask]))**2)


def fit_ccdf_params(theta_vals, ccdf_emp, theta_med):
    result = minimize(
        ccdf_loss,
        x0=[-np.log(2), 0.6],  # Initial guess for a and b
        bounds=[(-5, -0.001), (0.1, 2)],  # Reasonable bounds
        args=(theta_vals, ccdf_emp, theta_med),
        method='L-BFGS-B'
    )
    a_fit, b_fit = result.x
    print(f"Fitted CCDF parameters: a = {a_fit:.4f}, b = {b_fit:.4f}")
    return a_fit, b_fit

def fit_k_m(S_all, sizes_all, bin_centers, observed_medians):
    print('Total number of resolved sources:', len(sizes_all))

    result = minimize(
        loss,
        x0=[2, 0.3],  # Initial guess: k=2, m=0.3
        args=(bin_centers, observed_medians),
        bounds=[(0.1, 50), (0.01, 1)],
        method='L-BFGS-B' )

    k_fit, m_fit = result.x
    print(f"Fitted k: {k_fit:.3f}, m: {m_fit:.3f}")
    return k_fit, m_fit

def size_dist(S_all, sizes_all,flux_0,sizes_0, bin_centers, observed_medians, k_fit, m_fit):
    S_grid = np.logspace(np.log10(np.min(S_all)), np.log10(np.max(S_all)), 1000)
    theta_med_1 = median_size1(S_grid, k, m)
    theta_med_2 = median_size2(S_grid, k_fit, m_fit)
    plt.figure(figsize=(8, 6))
    plt.scatter(S_all, sizes_all, s=7, alpha=0.3, label='Resolved sources')
    plt.scatter(flux_0, sizes_0, s=15, alpha=0.3, marker='x',label='Unresolved sources',color='blue')
    plt.scatter(bin_centers, observed_medians, color='red', marker='*', s=100, label='Observed Medians')
    plt.plot(S_grid, theta_med_1, color='black', linewidth=2, label=r"$\Theta_{\mathrm{med,1}}$")
    plt.plot(S_grid, theta_med_2, color='purple', linewidth=2, label=r"$\Theta_{\mathrm{med,2}}$",)
    plt.plot(S_grid,theta_max(S_grid,sigma_low, bmaj, bmin),label="Maximum size", color="green",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_low, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]),label="Minimum size", color="orange",linestyle="--")
    plt.plot(S_grid,theta_max(S_grid,sigma_high, bmaj, bmin), color="green",linestyle="--")
    plt.plot(S_grid,theta_min(S_grid,sigma_high, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]), color="orange",linestyle="--")
    print('Min deconvolved size limit',np.min(theta_min(S_grid,sigma_high, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1])))
    plt.xscale('log')
    plt.ylim(-5,50)
    # plt.yscale('log')
    plt.xlabel(r'$S_T$ (mJy)', size=22)
    plt.ylabel(r'$\Theta$ (arcsec)', size=22)
    plt.legend(fontsize=15)
    #plt.grid(True, which='both', ls='--')
    plt.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
    plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/size_distribution.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()

def CCDF_plot(theta_med_1 ,theta_med_2,theta_vals, ccdf_emp,theta_lims_grid):
    ccdf_model_1 = integral_dist(theta_med_1, theta_lims_grid, a=a, b=b)
    ccdf_model_2 = integral_dist(theta_med_2, theta_lims_grid, a=a_fit, b=b_fit)

    theta_lims = theta_lim(S_grid, sigma_low, bmaj, bmin, alpha_beta_array[0], alpha_beta_array[1])
    ccdf_model_1_lim = integral_dist(theta_med_1, theta_lims, a=a_fit, b=b_fit)
    ccdf_model_2_lim = integral_dist(theta_med_2, theta_lims, a=a_fit, b=b_fit)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(theta_vals, ccdf_emp, drawstyle='steps-post', linewidth=3)
    ax1.plot(theta_lims_grid, ccdf_model_1, '-', label=r"$\Theta_{\mathrm{med,1}}$", linewidth=3,color='black')
    ax1.plot(theta_lims_grid, ccdf_model_2, '-', label=r"$\Theta_{\mathrm{med,2}}$", linewidth=3,color='purple')
    ax1.axvline(bmaj, ls='--', linewidth=2, color='black')
    ax1.set_xlabel(r"$\Theta_{\mathrm{lim}} \, \mathrm{[arcsec]}$", size=18)
    ax1.set_ylabel(r"$h(>\Theta_{\mathrm{lim}})$", size=18)
    ax1.set_xlim(5, 25)
    ax1.legend(fontsize=16)
    ax1.tick_params(axis='both', which='major', direction='in', length=4, width=1.5)
    ax1.tick_params(axis='both', which='minor', direction='in', length=8, width=1)

    c1 = 1 / (1 - ccdf_model_1_lim)
    c2 = 1 / (1 - ccdf_model_2_lim)
    ax2.plot(S_grid, c1, label=r"$\Theta_{\mathrm{med,1}}$", linewidth=3,color='black')
    ax2.plot(S_grid, c2, label=r"$\Theta_{\mathrm{med,2}}$", linewidth=3,color='purple')
    ax2.set_xlabel(r"$S_T$ [mJy]", size=18)
    ax2.set_ylabel(r"$c = 1/[1 - h(>\Theta_{\mathrm{lim}})]$", size=18)
    ax2.set_xscale('log')
    ax2.legend(fontsize=16)
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1.5)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolution_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.tight_layout()
    plt.show()

    interp_Func1 = interp1d(S_grid, c2, bounds_error=False, fill_value='extrapolate')
    with open('/home/kincaid/Desktop/Saraswati_codes/resolution_interp_func.pkl', 'wb') as f:
        pickle.dump(interp_Func1, f)
    return interp_Func1

    
if __name__=="__main__":
    cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/catalogs/MeerKAT_eddington_combined.fits')
    mask = cat['Maj'] == 0
    flux = cat['Total_flux'][~mask] * 1e3  # mJy
    breakpoint()
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
    #alpha_beta_array= [1.1, 1.35]
    alpha_beta_array= [1.1, 1.7] #[1.1, 1.35] #
    bin_centers1, medians1 = equal_source_bins(flux[mask_low], sizes[mask_low], nbins=5)
    bin_centers2, medians2 = equal_source_bins(flux[~mask_low], sizes[~mask_low], nbins=6)
    bin_centers = np.concatenate([bin_centers1, bin_centers2])
    median_sizes = np.concatenate([medians1, medians2])

    mask1= flux >1
    sizes=sizes[mask1]
    flux=flux[mask1]
    k_fit, m_fit = fit_k_m(flux, sizes, bin_centers, median_sizes)
    print(f"Fitted parameters: k = {k_fit:.4f}, m = {m_fit:.4f}")

    bin_size=50
    min_flux, max_flux = 0.01, 100
    S_grid = np.logspace(np.log10(min_flux), np.log10(max_flux),  bin_size)
    theta_lims_grid = np.logspace(np.log10(0.01), np.log10(50), bin_size)
    k,m=8,0.03
    theta_med_1 = median_size1(S_grid, k, m)
    theta_med_2 = median_size2(S_grid, k_fit, m_fit)
    theta_vals, ccdf_emp = complementary_ecdf(sizes, bin_size=bin_size)
    a, b = fit_ccdf_params(theta_vals, ccdf_emp, theta_med_1)
    #a,b=-np.log(2),0.62
    a_fit, b_fit = fit_ccdf_params(theta_vals, ccdf_emp, theta_med_2)

    interp_func=CCDF_plot(theta_med_1 ,theta_med_2,theta_vals, ccdf_emp,theta_lims_grid)

    sizes = np.sqrt(major * minor)
    sizes_0 = np.sqrt(major_0 * minor_0)
    flux = cat['Total_flux'][~mask] * 1e3  # mJy
    size_dist(flux, sizes,flux_0,sizes_0,  bin_centers, median_sizes, k_fit, m_fit)
  
    # Output
    