import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pickle

def binning(flux, nbins):
    pers = np.linspace(0, 100, nbins+1)
    #Range_x = np.percentile(flux, pers)
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max()), num=nbins+1)
    return(Range_x)

def median_size2(S, k, m):
    S = np.asarray(S)
    t = np.where(S < 1, k,  k * S**m)
    return t

def median_size1(S, k, m):
    S = np.asarray(S)
    t = k * S**m
    return t

def func(a,b,S,N):
    f=a+b/(S/N)
    return(f)

def integral_dist(theta_med,theta_lim,a=-np.log(2),b=0.62):
    x=np.exp(a*(theta_lim/theta_med)**b)
    return(x)

def theta_max(S,sigma, bmaj, bmin):
    theta_max=np.sqrt(bmaj*bmin)*np.sqrt((S/(threshold*sigma)))
    return(theta_max)

def theta_min(S,N, bmaj, bmin,a,b):
    mean_beam=np.sqrt(bmaj*bmin)
    theta_min=mean_beam*np.sqrt(func(a,b,S,N)-1)
    return(theta_min)

# def theta_lim(S, sigma, bmaj, bmin, a,b):
#         x=theta_min(S,sigma, bmaj, bmin, a,b)
#         y=theta_max(S, sigma, bmaj, bmin)
#         theta_lim=np.maximum(x,y)
#         mask=np.isnan(theta_lim)
#         theta_lim=theta_lim[~mask]
#         S=S[~mask]
#         return(S,theta_lim)

def theta_lim(S, sigma, bmaj, bmin, a,b):
        x=theta_min(S, sigma, bmaj, bmin, a,b)
        y=theta_max(S, sigma, bmaj, bmin)
        theta_lim=np.maximum(x,y)
        return(theta_lim)

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

def complementary_ecdf(data, bin_size):
    counts, bin_edges = np.histogram(data, bins=bin_size, density=False)
    cumsum = np.cumsum(counts[::-1])[::-1]  # reverse cumulative sum
    ccdf = cumsum / cumsum[0]  # normalize to [0, 1]
    return bin_edges[:-1], ccdf

def plot(true_sizes,sigma_high, alpha_beta_array):
    print('Number of sources for integral distribution A2631',len(true_sizes))
    S=np.linspace(0.01,100,10000)
    x1, y1 = complementary_ecdf(true_sizes,bin_size=40)
    fig, axes = plt.subplots(ncols=2, figsize=(17, 7))  # Create a side-by-side layout
    ax1 = axes[0]    
    ax1.plot(x1, y1, drawstyle='steps-post', linewidth=2, label=name, color='blue')
    breakpoint()
    Theta_lim= theta_lim(S, sigma_high, bmaj, bmin, alpha_beta_array[0], alpha_beta_array[1])

    #Theta_lim=np.linspace(1, 100, num=10000)
    h1= integral_dist(median_size1(S, k=8, m=0.3),  Theta_lim)
    h2= integral_dist(median_size2(S, k=8, m=0.03),  Theta_lim)
    ax1.plot(Theta_lim,h1 ,color='black',label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    ax1.plot(Theta_lim,h2,color='purple',label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    ax1.set_xlim(5, 30)
    ax1.set_ylabel(r"$h(>\Theta_{\mathrm{lim}})$", size=24)
    ax1.set_xlabel(r"$\Theta_{\mathrm{lim}} \, \mathrm{[arcsec]}$", size=24)
    ax1.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    ax1.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax1.set_xscale('linear')
    ax1.legend(fontsize=22)
    ax2 = axes[1]
    
    c1= 1/(1 - integral_dist(median_size1(S, k=8, m=0.3), Theta_lim))
    c2= 1/(1 - integral_dist(median_size2(S, k=8, m=0.03), Theta_lim))

    interp_Func1 = interp1d(S, c2, bounds_error=False, fill_value=0)
    with open('/home/kincaid/Desktop/Saraswati_codes/resolution_interp_func.pkl', 'wb') as f:
        pickle.dump(interp_Func1, f)

    ax2.plot(S, c1,color='black',label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    ax2.plot(S, c2,color='purple',label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    # ax2.plot(S, 1 / (1 - integral_dist(median_size(S, k=4, m=0.03), theta_lim(S,SN, sigma_high[0], bmaj, bmin, alpha_beta_array[0][0], alpha_beta_array[0][1]))),color='black',label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    #ax2.plot(S, 1 / (1 - integral_dist(median_size(S, k=4, m=0.3), theta_lim(S,SN, sigma_high[1], bmaj, bmin, alpha_beta_array[0][0], alpha_beta_array[0][1]))),color='purple',label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    ax2.set_ylabel(r"c=$1/[1-h(>\Theta_{\mathrm{lim}})$]", size=24)
    ax2.set_xlabel(r"$S_T$ [mJy]", size=24)
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax2.set_xscale('log')
    ax2.legend(fontsize=22)
    plt.xlim(0.1,100)
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolution_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
    return interp_Func1

def correction(cat,interp_func):
    mask=cat['Maj'] !=0
    flux=cat['Total_flux'][mask]*1e3
    corr=interp_func(np.array(flux))
    S_corr=flux*1/corr
    cat['Total_flux'][mask] = S_corr*1e-3  # Convert back to Jy
    outname=f'{path}{name}/catalogs/{name}_resolution_corr_srl.fits'
    cat.write(outname, overwrite=True)
    print(f"Resolution bias corrected catalog written to {outname}")

if "__main__":
    names=['A2631','Zwcl2341']
    bmaj=8.7159409207893
    bmin=7.209874015964243
    threshold=5
    path='/home/kincaid/Desktop/Saraswati_codes/'
    true_sizes_array=[]
    theta_med_array=[]
    m_array=[]
    fluxes_array=[]
    alpha_beta_array=[[1.1,1.35],[1.1,1.7]]
    #alpha_beta_array=[1.1,1.7]
    sigma_high_array=[0.06,0.06]
    min_flux_array=[0.2,0.2]
    max_flux_array=[1000,1000]
    
    for i,name in enumerate(names):
        cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_eddington_corr_srl.fits')
        print('This is',i)
        true_sizes,fluxes=info(cat,min_flux_array[i],max_flux_array[i])
        true_sizes_array.append(true_sizes)
        fluxes_array.append(fluxes)
        interp_func=plot(true_sizes_array[i],sigma_high_array[i], alpha_beta_array[i])
        correction(cat,interp_func)
    

    