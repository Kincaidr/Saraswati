import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from astropy.io import fits
from astropy.table import Table
import json



# def func(a,b,x):
#     f=a/(1+b/x)
#     return(f)

# def func(a,b,x):
#     f=1+a*x**b
#     return(f)

def func(a,b,x):
    f=a+b/x
    return(f)

def constraint(params, x,y):
    a, b = params
    f=func(a,b,x)
    counter = np.sum( (y > f))
    val = (0.95 - (counter/len(x)))**2
    #breakpoint()
    return val

def optimization(catalog,rms,initial_cond, method):
   
    cat= Table.read(catalog)
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    rms=cat['Isl_rms']
    y=flux_total/flux_peak
    mask=y <= 1
    y=y[mask]
    x=flux_peak/rms
    x=x[mask]
    order=np.argsort(x)
    x, y = x[order], y[order]
    minimizer_kwargs = {"args":(x,y)}
    res = minimize(constraint, initial_cond, args=minimizer_kwargs["args"], method=method )
    x0=res.x[0]
    x1=res.x[1] 
    results=np.array([x0,x1])
    print('success',res.success)
    print('status',res.status)
    print('message',res.message)
    print('nit',res.nit)
    print('x0',x0)
    print('x1',x1)
    return(results)

def plot_curve(catalog,rms, results, method, output_plot):

    cat= Table.read(catalog)
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    rms=cat['Isl_rms']
    x0=results[0]
    x1=results[1] 
    y=flux_total/flux_peak
    x=flux_peak/rms
    source_tot=len(x)
    t = np.logspace(0.1, 3)
    curve_theory=func(x0,x1,t)
    curve_theory_inverse=1/curve_theory
    curve_real=func(x0,x1,x)
    curve_real_inverse=1/curve_real
    unresolved_mask= (y < curve_real) & (y >curve_real_inverse)
    resolved_mask= ~unresolved_mask
    total_unresolved=np.sum(unresolved_mask)
    total_resolved=np.sum(resolved_mask)
    cond=(y <= 1) & (y > curve_real_inverse)
    num=np.sum(cond)/np.sum(y <=1)
    print('Percentage',num*100)
    print('Method', method,'function')
    print("Total sources",source_tot)
    print("unresolved fraction",(total_unresolved/source_tot)*100)
    print("resolved fraction",(total_resolved/source_tot)*100)
    
    #fig,ax = plt.figure(figsize=(13, 10))
    fig,ax=plt.subplots(figsize=(8, 6))
    #ax.fill_between(t,curve_theory,curve_theory_inverse,alpha=0.6,label='unresolved sources')
    plt.plot(t,curve_theory,linewidth=2,color='black')
    plt.plot(t,curve_theory_inverse,linewidth=2,color='black')
    plt.axhline(y=1,color='r',linestyle='--')
    plt.scatter(x[unresolved_mask],y[unresolved_mask],alpha=0.5,s=12,color='green')
    plt.scatter(x[resolved_mask],y[resolved_mask],alpha=0.5,s=12,color='blue')
    plt.xlabel(r'$S_P/\sigma$',fontsize=20)
    plt.ylabel(r'$S_T/ S_P$',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=18, length=5, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=18, length=5, width=1)
    plt.xscale('log')
    plt.yscale('log')
    #plt.title( rf'$solver = {method}, resolved = {total_resolved}, unresolved = {total_unresolved}, a = {x0:.3f}, b = {x1:.3f}$', fontsize=10) 
    #plt.title(rf'$solver = {method}, resol = {total_resolved}, unresol = {total_unresolved}, a = {x0:.3f}, b = {x1:.3f}$')
    plt.xlim([3, 1000])
    plt.ylim([0.3, 100])
    plt.tight_layout()
    plt.legend()
    output_name=name+'_resolved_unresolved.png'
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/'+output_name)
    print(output_plot+output_name + ' saved')


def flux_correction(catalog,rms, results, output_cat):

    cat= Table.read(catalog)
    x0=results[0]
    x1=results[1] 
    flux_total=cat['Total_flux']
    flux_peak=cat['Peak_flux']
    y=flux_total/flux_peak
    x=flux_peak/rms
    curve_real=func(x0,x1,x)
    curve_real_inverse=1/curve_real
    unresolved_mask= (y < curve_real) & (y >curve_real_inverse)
    cat['Maj'][unresolved_mask]=0
    cat['Min'][unresolved_mask]=0
    cat['Total_flux'][unresolved_mask]=cat['Peak_flux'][unresolved_mask]
    corrected_cat=name+'_srl_flux_corr.fits'
    cat.write(output_cat+corrected_cat,format='fits', overwrite=True)
    
def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

if __name__ == "__main__":

    name='Zwcl2341'
    catalogs="/home/kincaid/Desktop/Saraswati_codes/"+name+"/catalogs/"
    plots="/home/kincaid/Desktop/Saraswati_codes/"+name+"/plots/"
    rms = 15e-6
    method='Powell'#Powell','CG','BFGS','L-BFGS-B','TNC','COBYLA', 'Nelder-Mead'] 
    #initial_cond=[3.2,-0.9]
    #initial_cond=[1.09,2.7]  #99%
    #initial_cond=[1.05,1.35]  #95% A2631
    initial_cond=[1,1.7] #95% zwcl
    output_cat=catalogs
    output_plot=plots
    real_catalog=catalogs+name+"_full_srl.fits"
    catalogs=[real_catalog]
    for catalog in catalogs:
        results = optimization(catalog,rms,initial_cond, method)
        plot_curve(catalog,rms, results, method, output_plot)
        flux_correction(catalog,rms, results, output_cat)




