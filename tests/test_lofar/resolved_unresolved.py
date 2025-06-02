
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from astropy.io import fits
from astropy.table import Table
import json
import os

# def func(a,b,x):
#     f=a/(1+b/x)
#     return(f)

def func(a,b,x):
    f=1+a*x**b
    return(f)

# def func(a,b,x):
#     f=a+b/x
#     return(f)

def constraint(params, x,y):
    a, b = params
    f=func(a,b,x)
    counter = np.sum( (y > f))
    val = abs(0.95 - (counter/len(x)))
    return val

def optimization(catalog,rms,initial_cond, method):
   
    cat= Table.read(catalog)
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    y=flux_total/flux_peak
    mask=y <= 1
    y=y[mask]
    x=flux_peak/rms
    x=x[mask]
    order=np.argsort(x)
    x, y = x[order], y[order]
    minimizer_kwargs = {"args":(x,y)}
    res = minimize(constraint, initial_cond, args=minimizer_kwargs["args"], method=method)
    x0=res.x[0]
    x1=res.x[1] 
    results=np.array([x0,x1])
    print('success',res.success)
    print('status',res.status)
    print('message',res.message)
    print('nit',res.nit)
    return(results)

def plot_curve(catalog,rms, results, method, output_img):

    cat= Table.read(catalog)
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    x0=results[0]
    x1=results[1] 
    y=flux_total/flux_peak
    x=flux_peak/rms
    source_tot=len(x)
    breakpoint()
    t = np.logspace(0.1, 3)
    curve=func(x0,x1,t)
    curve_inverse=1/curve
    curve1=func(x0,x1,x)
    curve1_inverse=1/curve1
    unresolved_mask= (y > curve1) & (y <curve1_inverse)
    resolved_mask= ~unresolved_mask
    total_unresolved=np.sum(unresolved_mask)
    total_resolved=np.sum(resolved_mask)
    cond=(y <= 1) & (y > curve1)
    num=np.sum(cond)/np.sum(y <=1)
    print('Percentage',num*100)
    print('Method', method,'function')
    print("Total sources",source_tot)
    print("Total unresolved sources",total_unresolved)
    print("Total resolved sources",total_resolved)
    fig,ax=plt.subplots()
    ax.fill_between(t,curve,curve_inverse,alpha=0.6,label='unresolved sources')

    plt.plot(t,curve,linewidth=2,color='black')
    plt.plot(t,curve_inverse,linewidth=2,color='black')
    plt.axhline(y=1,color='r',linestyle='--')
    plt.scatter(x,y,alpha=0.5,s=10,color='grey')
    plt.xlabel(r'$S_P/\sigma$',fontsize=18)
    plt.ylabel(r'$S_T/ S_P$',fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.title( rf'$solver = {method}, resolved = {total_resolved}, unresolved = {total_unresolved}, a = {x0:.3f}, b = {x1:.3f}$', fontsize=10) 
    #plt.title(rf'$solver = {method}, resol = {total_resolved}, unresol = {total_unresolved}, a = {x0:.3f}, b = {x1:.3f}$')
    # plt.xlim([5, 1000])
    # plt.ylim([-10, 10])
    plt.tight_layout()
    plt.legend()
    _,save_name=os.path.split(catalog)
    output=output_img+save_name.replace('.fits','_resolved_unresolved.png')
    plt.savefig(output)
    print(output + ' saved')


def flux_correction(catalog,rms, results, output_cat):

    cat= Table.read(catalog)
    x0=results[0]
    x1=results[1] 
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    y=flux_total/flux_peak
    x=flux_peak/rms
    curve1=func(x0,x1,x)
    curve1_inverse=1/curve1
    unresolved_mask= (y > curve1) & (y <curve1_inverse)
    cat['Total_flux'][unresolved_mask]=cat['Peak_flux'][unresolved_mask]
    corrected_cat='merged_rec_flux_corr.fits'
    cat.write(output_cat+corrected_cat,format='fits', overwrite=True)
    print(output_cat + ' written')
    
def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

if __name__ == "__main__":

    config_path = 'A2631.json'
    config = read_config(config_path)
    path=config['path']
    catalogs=config['catalogs']
    plots=path+config['plots']
    rec_catalog=path+catalogs+config['merged_rec']
    rms = config['sigma']*10**-3
    method='Nelder-Mead'#Powell','CG','BFGS','L-BFGS-B','TNC','COBYLA', 'Nelder-Mead'] 
    initial_cond=[3.2,-0.9]
    

    results = optimization(rec_catalog,rms,initial_cond, method)
    plot_curve(rec_catalog,rms, results, method, plots)
    flux_correction(rec_catalog,rms, results, catalogs)
