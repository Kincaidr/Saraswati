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


def info(cat, min_flux, max_flux):
    mask=(cat['Maj']==0)
    flux=cat['Total_flux'][~mask]*1e3
    min=cat['Min'][~mask]*3600
    maj=cat['Maj'][~mask]*3600
    size=np.sqrt(maj*min)
    nbins=20
    Range_x = binning(flux, nbins)
    bin_centers=[]
    median_sizes=[]

    for i in range(len(Range_x)-1):
        bin_mask = (flux >= Range_x[i]) & (flux < Range_x[i+1])
        if np.any(bin_mask):
            bin_centers.append((Range_x[i] + Range_x[i+1]) / 2)
            median_sizes.append(np.median(size[bin_mask]))

    bin_centers = np.array(bin_centers)
    median_sizes = np.array(median_sizes)
    popt, pcov = curve_fit(median_size, bin_centers, median_sizes, p0=[10, 0.1])
    mask= (min_flux <flux) & (flux < max_flux)
    fluxes=flux[mask]
    true_sizes=size[mask]
    theta_med=popt[0]
    m=popt[1]
    return(true_sizes,fluxes,theta_med,m,thresh)


def plot(true_sizes,theta_med,m,sigma_high, alpha_beta_array):
    #S_new=((thetas/np.sqrt(bmaj*bmin))**2+1)*5*sigma_high
    print('Number of sources for integral dsitribution A2631',len(true_sizes[0]))
    print('Number of sources for integral dsitribution Zwcl2341',len(true_sizes[1]))
    S=np.linspace(1,100,num=1000)
    thetas=np.linspace(0,120, num=1000)
    fig, axes = plt.subplots(ncols=2, figsize=(20, 8))  # Create a side-by-side layout
    ax1 = axes[0]
    ax1.ecdf(true_sizes[0], complementary=True )
    ax1.ecdf(true_sizes[1], complementary=True )
    ax1.plot(thetas, integral_dist(median_size(S, theta_med[0], m[0]), thetas), label='A2631',color='blue',linestyle='--')
    ax1.plot(thetas, integral_dist(median_size(S, theta_med[1], m[1]), thetas), label='Zwcl2341',color='red',linestyle='--')
    #ax1.plot(thetas, integral_dist(median_size(S_new, popt[0], popt[1]), thetas), label='Median sizes new')
    #ax1.plot(thetas, integral_dist(theta_med, thetas), label='Median size fixed')
    ax1.set_xlim(0, 100)
    ax1.set_ylabel(r"$h(>\Theta_{\mathrm{lim}})$", size=22)
    ax1.set_xlabel(r"$\Theta_{\mathrm{lim}} \, \mathrm{[arcsec]}$", size=22)
    ax1.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=16)
    ax1.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax1.set_xscale('linear')
    ax1.minorticks_on()
    ax1.legend(fontsize=18)

    ax2 = axes[1]
    S=np.linspace(0.001,1000,num=1000)
    #ax2.plot(S_new, 1 / (1 - integral_dist(median_size(S_new, popt[0], popt[1]), thetas)))
    ax2.plot(S, 1 / (1 - integral_dist(median_size(S, theta_med[0], m[0]), theta_lim(S, sigma_high[0], bmaj, bmin, alpha_beta_array[0][0], alpha_beta_array[0][1]))),color='blue',label='A2631')
    ax2.plot(S, 1 / (1 - integral_dist(median_size(S, theta_med[1], m[1]), theta_lim(S, sigma_high[1], bmaj, bmin, alpha_beta_array[1][0], alpha_beta_array[1][1]))),color='red',label='Zwcl2341')
    ax2.set_ylabel(r"c=$1/[1-h(>\Theta_{\mathrm{lim}})$]", size=22)
    ax2.set_xlabel(r"$S_T$ [mJy]", size=22)
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=16)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=16)
    ax2.set_xscale('log')
    ax2.minorticks_on()
    ax2.legend(fontsize=18)
    plt.xlim(0.001,1000)
    # Adjust layout for better spacing
    plt.tight_layout()
    plt.savefig('plots/resolution_bias.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()

def correction(cat,fluxes,theta_med,sigma_high, alpha_betas):
    S=np.linspace(0.001,1000,num=1000)
    corr=1/(1 - integral_dist(median_size(S, theta_med, m), theta_lim(S, sigma_high, bmaj, bmin, alpha_betas[0],alpha_betas[1])))
    output_file = name+"_resolution_bias_correction.txt"
    with open(output_file, "w") as file:
        for x ,y in zip(corr, S):
            file.write(f"{x} {y} \n")

if "__main__":
    names=['A2631','Zwcl2341']
    bmaj=8.7159409207893
    bmin=7.209874015964243
    thresh=5
    plots='/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'
    true_sizes_array=[]
    theta_med_array=[]
    m_array=[]
    fluxes_array=[]
    alpha_beta_array=[[1.1,1.35],[1.1,1.7]]
    sigma_high_array=[0.04,0.04]
    min_flux_array=[1,1]
    max_flux_array=[100,100]
    for i,name in enumerate(names):
        cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_srl_flux_corr.fits')
        print('This is',i)
        true_sizes,fluxes,theta_med,m,thresh=info(cat,min_flux_array[i],max_flux_array[i])
        true_sizes_array.append(true_sizes)
        theta_med_array.append(theta_med)
        m_array.append(m)
        fluxes_array.append(fluxes)
        correction(cat,fluxes,theta_med,sigma_high_array[i], alpha_beta_array[i])
        print(alpha_beta_array[i][0],alpha_beta_array[i][1])
    plot(true_sizes_array,theta_med_array,m_array,sigma_high_array, alpha_beta_array)

    