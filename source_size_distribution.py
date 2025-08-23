import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(a,b,S,N):
    f=a+b/(S/N)
    return(f)

def integral_dist(theta_med,theta_lim,a=-np.log(2),b=0.62):
    x=np.exp(b*(theta_lim/theta_med)**a)
    return(x)

def theta_max(S,N, bmaj, bmin):
    theta_max=np.sqrt(bmaj*bmin)*np.sqrt((S/(threshold*N)))
    return(theta_max)

def theta_min(S,N, bmaj, bmin,a,b):
    mean_beam=np.sqrt(bmaj*bmin)
    theta_min=mean_beam*np.sqrt(func(a,b,S,N)-1)
    return(theta_min)

def median_size2(S, k, m):
    S = np.asarray(S)
    t = np.where(S < 1, k,  k * S**m)
    return t

def median_size1(S, k, m):
    S = np.asarray(S)
    t = k * S**m
    return t

def equal_source_bins(flux,size,nbins):
    bin_centers = []
    median_sizes = []
    percentiles = np.linspace(0, 100, nbins + 1)
    bin_edges = np.percentile(flux, percentiles)

    for i in range(len(bin_edges) - 1):
        bin_mask = (flux >= bin_edges[i]) & (flux < bin_edges[i+1])
        if i == len(bin_edges) - 2:
            bin_mask = (flux >= bin_edges[i]) & (flux <= bin_edges[i+1])

        if np.any(bin_mask):
            bin_centers.append(np.median(flux[bin_mask]))
            median_sizes.append(np.median(size[bin_mask]))
    return np.array(bin_centers), np.array(median_sizes)


def plot(cat, bmaj, bmin, sigma_high,sigma_low,axs,color,alpha_beta_array):
    total_sources=len(cat['Maj'])
    print('Total number of sources', total_sources)
    mask=(cat['Maj']==0)
    flux=cat['Total_flux'][~mask]*1e3
    min=cat['Min'][~mask]*3600
    maj=cat['Maj'][~mask]*3600
    size=np.sqrt(maj*min)

    mask1=flux <1
    bin_centers1,median_sizes1= equal_source_bins(flux[mask1],size[mask1],nbins=5)
    bin_centers2,median_sizes2= equal_source_bins(flux[~mask1],size[~mask1],nbins=6)

    bin_centers = np.concatenate((bin_centers1, bin_centers2))
    median_sizes = np.concatenate((median_sizes1, median_sizes2))
    #popt, pcov = curve_fit(median_size, bin_centers, median_sizes, p0=[10, 0.1])
    #k_fit, m_fit = popt
    #print(f"Fitted parameters: k = {k_fit}, m = {m_fit}")
    S=np.logspace(-1.5,2.5,10000)
    flux=cat['Total_flux']*1e3
    min=cat['Min']*3600
    maj=cat['Maj']*3600
    size=np.sqrt(maj*min)
    axs.scatter(flux, size, s=10, alpha=0.3,color=color)
    axs.scatter(bin_centers, median_sizes, label='Median size', marker='*',color=color_median,s=100)
    axs.plot(S,theta_max(S,sigma_low, bmaj, bmin),label="Maximum size", color="green",linestyle="--")
    axs.plot(S,theta_min(S,sigma_low, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]),label="Minimum size", color="orange",linestyle="--")

    axs.plot(S,theta_max(S,sigma_high, bmaj, bmin), color="green",linestyle="--")
    axs.plot(S,theta_min(S,sigma_high, bmaj, bmin,alpha_beta_array[0],alpha_beta_array[1]), color="orange",linestyle="--")
    axs.plot(S, median_size1(S, k=8, m=0.3), label=r"$\Theta_{\mathrm{med,1}}$", color="black",linewidth=2)
    axs.plot(S, median_size2(S, k=0.532, m=0.476), label=r"$\Theta_{\mathrm{med,2}}$", color="purple",linewidth=2)
    axs.set_xscale('log')
    axs.set_ylim(-2,80)
    axs.set_xlabel(r"$S_T$ [mJy]", size=22)
    axs.set_ylabel(r"$\Theta$ [arcsec]", size=22)
    axs.tick_params(axis='both', which='both', direction='in', length=6 , labelsize=14)  # 'in' means they point inward
    axs.tick_params(which='minor', length=3)
    #plt.title("Disitrubtion of real soure sizes with fitted median size model")


if "__main__":
    names=['A2631','Zwcl2341']
    bmaj=8.7159409207893
    bmin=7.209874015964243
    sigma_low=[0.016, 0.011]
    sigma_high=[0.06, 0.06]
    alpha_beta_array= [[1.1, 1.35], [1.1, 1.7]]

    #alpha_beta_array=[1.1,1.7] #Zwcl2341
    colors=['blue','red']
    color_median='cyan'
    threshold=5
    fig, axs = plt.subplots(1, 2, figsize=(8, 6), sharey=True)
    for i,name in enumerate(names):
        cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_flux_corr_srl.fits')
        plot(cat, bmaj, bmin, sigma_high[i],sigma_low[i],axs[i],colors[i],alpha_beta_array[i])
    plt.legend(fontsize=15)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/size_distribution.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()