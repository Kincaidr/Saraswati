import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(a,b,S,N,t):
    f=a+b/(S/(t*N))
    return(f)

def integral_dist(theta_med,theta_lim,a=-np.log(2),b=0.62):
    x=np.exp(b*(theta_lim/theta_med)**a)
    return(x)

def theta_max(S, sigma, bmaj, bmin):
    theta_max=np.sqrt(bmaj*bmin)*np.sqrt((S/(threshold*sigma))-1)
    return(theta_max)

def theta_min(S,sigma, bmaj, bmin,alpha_betas):
    t=5
    theta_min=np.sqrt(bmaj*bmin)*np.sqrt(func(alpha_betas[0],alpha_betas[1],S,sigma,t)-1)
    return(theta_min)

# def theta_min(S,sigma, bmaj, bmin,a=1.04,b=2.6):
#     theta_min=np.sqrt(bmaj*bmin)*np.sqrt(1/(a+b/(S/4*sigma))-1)
#     return(theta_min)

def median_size( S,k,m):
    t=k*S**m
    return(t)

def binning(flux, nbins):
    pers = np.linspace(0, 100, nbins+1)
    #Range_x = np.percentile(flux, pers)
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max()), num=nbins+1)
    return(Range_x)


def plot(cat, bmaj, bmin, sigma_low, sigma_high):
    total_sources=len(cat['Maj'])
    print('Total number of sources', total_sources)
    mask=(cat['Maj']==0)
    flux=cat['Total_flux'][~mask]*1e3
    min=cat['Min'][~mask]*3600
    maj=cat['Maj'][~mask]*3600
    size=np.sqrt(maj*min)
    nbins=10
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
    k_fit, m_fit = popt
    print(f"Fitted parameters: k = {k_fit}, m = {m_fit}")
    S=np.logspace(-1.5,2,10000)
    flux=cat['Total_flux']*1e3
    min=cat['Min']*3600
    maj=cat['Maj']*3600
    size=np.sqrt(maj*min)

    fig, ax1 = plt.subplots(figsize=(8, 8))
    plt.scatter(flux, size, s=10, alpha=0.3,color=color)
    plt.scatter(bin_centers, median_sizes, label='Median size', marker='*',color=color_median)
    plt.plot(S,theta_max(S,sigma_low, bmaj, bmin),label="Maximum size", color="green",linestyle="--")
    plt.plot(S,theta_min(S,sigma_low, bmaj, bmin,alpha_beta_array),label="", color="orange",linestyle="--")
    plt.plot(S,theta_max(S,sigma_high, bmaj, bmin),label="", color="green",linestyle="--")
    plt.plot(S,theta_min(S,sigma_high, bmaj, bmin,alpha_beta_array),label="Minimum size", color="orange",linestyle="--")
    plt.plot(S, median_size(S, *popt), label="Windhorst fitted Curve", color="purple")
    plt.xscale('log')
    plt.ylim(-2,80)
    plt.xlabel(r"$S_T$ [mJy]", size=22)
    plt.ylabel(r"$\Theta$ [arcsec]", size=22)
    ax1.tick_params(axis='both', which='both', direction='in', length=6 , labelsize=14)  # 'in' means they point inward
    ax1.tick_params(which='minor', length=3)
    #plt.title("Disitrubtion of real soure sizes with fitted median size model")
    plt.legend(fontsize=16)
    plt.savefig('plots/'+name+'_size_distribution.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()

if "__main__":
    name='A2631'
    cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_srl_flux_corr.fits')
    outname='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    bmaj=8.7159409207893
    bmin=7.209874015964243
    sigma_low=0.03
    sigma_high=0.3
    #alpha_beta_array=[1.05,1.045] #A2631
    alpha_beta_array=[1,1.7] #Zwcl2341
    color='blue'
    color_median='cyan'
    threshold=5
    plot(cat, bmaj, bmin, sigma_low, sigma_high)