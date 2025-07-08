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

# def constraint(params, x, y):
#     a, b = params
#     f = func(a, b, x)
#     diff = y - f
#     penalty = np.mean(np.maximum(diff, 0))  # penalize values above the curve
#     target = 0.95
#     coverage = np.sum(y <= f) / len(x)
#     print(f"Trying a={a:.4f}, b={b:.4f} => coverage={np.sum(y <= func(a,b,x)) / len(x):.4f}")
#     return (coverage - target)**2 + 0.1 * penalty

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

def plot_curve(catalog, rms, results, method, ax, color, label):
    from astropy.table import Table
    import numpy as np

    cat = Table.read(catalog)
    flux_total = cat['Total_flux']
    flux_peak = cat['Peak_flux']
    rms = cat['Isl_rms']
    
    x0, x1 = results
    y = flux_total / flux_peak
    x = flux_peak / rms
    source_tot = len(x)

    t = np.logspace(0.1, 3)
    curve_theory = func(x0, x1, t)
    curve_theory_inverse = 1 / curve_theory
    curve_real = func(x0, x1, x)
    curve_real_inverse = 1 / curve_real

    unresolved_mask = (y < curve_real) & (y > curve_real_inverse)
    resolved_mask = ~unresolved_mask

    total_unresolved = np.sum(unresolved_mask)
    total_resolved = np.sum(resolved_mask)

    cond = (y <= 1) & (y > curve_real_inverse)
    num = np.sum(cond) / np.sum(y <= 1)

    # Print stats
    print('Percentage', num * 100)
    print('Method', method, 'function')
    print("Total sources", source_tot)
    print("Unresolved fraction", (total_unresolved / source_tot) * 100)
    print("Resolved fraction", (total_resolved / source_tot) * 100)

    # Plot using passed axis
    ax.plot(t, curve_theory, linewidth=2, color='black')
    ax.plot(t, curve_theory_inverse, linewidth=2, color='black')
    ax.axhline(y=1, color='red', linestyle='--')

    ax.scatter(x[unresolved_mask], y[unresolved_mask], alpha=0.5, s=12, color='green')
    ax.scatter(x[resolved_mask], y[resolved_mask], alpha=0.5, s=12, color=color)
    ax.set_title(label, fontsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([3, 1000])
    ax.set_ylim([0.3, 100])
    ax.set_xlabel(r'$S_P/\sigma$', fontsize=22)
    ax.set_ylabel(r'$S_T/ S_P$', fontsize=22)
    ax.tick_params(axis='both', which='major', labelsize=12, length=5, width=1)
    ax.tick_params(axis='both', which='minor', labelsize=12, length=5, width=1)
    ax.legend(fontsize=10)


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
    # name='A2631'
    # catalogs="/home/kincaid/Desktop/Saraswati_codes/"+name+"/catalogs/"
    # plots="/home/kincaid/Desktop/Saraswati_codes/"+name+"/plots/"

    # rms = 16e-6
    # color='blue'
    # method='Powell'#Powell','CG','BFGS','L-BFGS-B','TNC','COBYLA', 'Nelder-Mead'] 
    # initial_cond=[1.03,1.4]  #95% A2631
    # #initial_cond=[1.01,1.7] #95% zwcl
    # output_cat=catalogs
    # output_plot=plots
    # real_catalog=catalogs+name+"_cut_srl.fits"
    # catalogs=[real_catalog]
    # for catalog in catalogs:
    #     results = optimization(catalog,rms,initial_cond, method)
    #     plot_curve(catalog,rms, results, method, output_plot)
    #     flux_correction(catalog,rms, results, output_cat)
    names = ['A2631', 'Zwcl2341']
    colors = ['blue', 'red']
    initial_conds = [[1.03, 1.4], [1.01, 1.7]]
    rms_values = [16e-6, 11e-6]  # Adjust if needed

    fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharey=True)

    for i, name in enumerate(names):
        catalogs_path = f"/home/kincaid/Desktop/Saraswati_codes/{name}/catalogs/"
        catalog = catalogs_path + name + "_cut_srl.fits"
        results = optimization(catalog, rms_values[i], initial_conds[i], method='Powell')
        plot_curve(catalog, rms_values[i], results, method='Powell', ax=axs[i], color=colors[i], label=names[i])
        flux_correction(catalog, rms_values[i], results, catalogs_path)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolved_unresolved.png', dpi=300)
    print('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/resolved_urnresolved.png')
    plt.show()




