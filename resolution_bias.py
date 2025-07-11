import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def binning(flux, nbins):
    pers = np.linspace(0, 100, nbins+1)
    #Range_x = np.percentile(flux, pers)
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max()), num=nbins+1)
    return(Range_x)

def median_size(S,k,m):
    t=k*S**m
    return(t)

def func(a,b,S,N,t):
    f=a+b/(S/(t*N))
    return(f)

def integral_dist(theta_med,theta_lim,a=-np.log(2),b=0.62):
    x=np.exp(a*(theta_lim/theta_med)**b)
    return(x)

def theta_max(S, sigma, bmaj, bmin):
    t=thresh
    theta_max=np.sqrt(bmaj*bmin)*np.sqrt((S/(t*sigma))-1)
    return(theta_max)

def theta_min(S,sigma, bmaj, bmin,a,b):
    t=thresh
    theta_min=np.sqrt(bmaj*bmin)*np.sqrt(func(a,b,S,sigma,t)-1)
    return(theta_min)

def theta_lim(S, sigma, bmaj, bmin, a,b):
    theta_lim_arr=[]
    for i in range(len(S)):
        x=theta_min(S[i],sigma, bmaj, bmin, a,b)
        y=theta_max(S[i], sigma, bmaj, bmin)
        theta_lim=max(x,y)
        theta_lim_arr.append(theta_lim)
    return(theta_lim_arr)

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


def plot(true_sizes,sigma_high, alpha_beta_array):
    #S_new=((thetas/np.sqrt(bmaj*bmin))**2+1)*5*sigma_high
    print('Number of sources for integral dsitribution A2631',len(true_sizes[0]))
    print('Number of sources for integral dsitribution Zwcl2341',len(true_sizes[1]))
    S=np.linspace(0.01,100,num=100000)
    theta_lim_arr=np.linspace(0,120, num=100000)
    fig, axes = plt.subplots(ncols=2, figsize=(20, 8))  # Create a side-by-side layout
    ax1 = axes[0]
    ax1.ecdf(true_sizes[0], complementary=True,linewidth=2,label='A2631',color='blue')
    ax1.ecdf(true_sizes[1], complementary=True,linewidth=2,label='Zwcl2341',color='red') 
    h1= integral_dist(median_size(S, k=4, m=0.03), theta_lim_arr)
    h2= integral_dist(median_size(S, k=4, m=0.0), theta_lim_arr)
    ax1.plot(S,h1 ,color='black',label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    ax1.plot(S,h2,color='purple',label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    ax1.set_xlim(0, 40)
    ax1.set_ylabel(r"$h(>\Theta_{\mathrm{lim}})$", size=22)
    ax1.set_xlabel(r"$\Theta_{\mathrm{lim}} \, \mathrm{[arcsec]}$", size=22)
    ax1.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    ax1.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax1.set_xscale('linear')
    ax1.minorticks_on()
    ax1.legend(fontsize=22)

    S=np.linspace(0.01,1000,num=100000)
    ax2 = axes[1]    
    breakpoint()
    ax2.plot(S, 1 / (1 - integral_dist(median_size(S, k=4, m=0.03), theta_lim(S, sigma_high[0], bmaj, bmin, alpha_beta_array[0][0], alpha_beta_array[0][1]))),color='black',label=r"$\Theta_{\mathrm{med,1}}$",linewidth=3)
    ax2.plot(S, 1 / (1 - integral_dist(median_size(S, k=4, m=0.3), theta_lim(S, sigma_high[1], bmaj, bmin, alpha_beta_array[1][0], alpha_beta_array[1][1]))),color='purple',label=r"$\Theta_{\mathrm{med,2}}$",linewidth=3)
    
    ax2.set_ylabel(r"c=$1/[1-h(>\Theta_{\mathrm{lim}})$]", size=22)
    ax2.set_xlabel(r"$S_T$ [mJy]", size=22)
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax2.set_xscale('log')
    ax2.minorticks_on()
    ax2.legend(fontsize=22)
    plt.xlim(0.01,1000)
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolution_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()


def correction(cat,sigma_high, alpha_betas):
    mask=cat['Maj'] !=0
    flux=cat['Total_flux'][mask]*1e3
    corr=1/(1 - integral_dist(median_size(flux, k=4, m=0.03), theta_lim(flux, sigma_high, bmaj, bmin, alpha_betas[0],alpha_betas[1])))
    S_corr=flux*corr
    cat['Total_flux'][mask] = S_corr
    outname=f'{path}{name}/catalogs/{name}_resolution_corr_srl.fits'
    cat.write(outname, overwrite=True)
    print(f"Resolution bias corrected catalog written to {outname}")

if "__main__":
    names=['A2631','Zwcl2341']
    bmaj=8.7159409207893
    bmin=7.209874015964243
    thresh=5
    path='/home/kincaid/Desktop/Saraswati_codes/'
    true_sizes_array=[]
    theta_med_array=[]
    m_array=[]
    fluxes_array=[]
    alpha_beta_array=[[1.1,1.35],[1.1,1.7]]
    sigma_high_array=[0.04,0.04]
    min_flux_array=[0,0]
    max_flux_array=[1000,1000]
    for i,name in enumerate(names):
        cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_eddington_corr_srl.fits')
        print('This is',i)
        true_sizes,fluxes=info(cat,min_flux_array[i],max_flux_array[i])
        true_sizes_array.append(true_sizes)
        fluxes_array.append(fluxes)
        #correction(cat,sigma_high_array[i], alpha_beta_array[i])
        print(alpha_beta_array[i][0],alpha_beta_array[i][1])
    plot(true_sizes_array,sigma_high_array, alpha_beta_array)

    