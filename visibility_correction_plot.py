import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def get_data(file1):
    table1=np.loadtxt(file1)
    bin1=table1[:,0]
    comp1=table1[:,1]
    comp1_corr=table1[:,2]
    err1=table1[:,3]
    return(bin1,comp1,comp1_corr,err1)

def plot(bins,comps,comp_corrs,errs):
    plt.errorbar(bins[0],comps[0],yerr=errs[0], label='A2631', color='blue')
    plt.errorbar(bins[0],comp_corrs[0],yerr=errs[0], color='blue',linestyle='--')
    plt.errorbar(bins[1],comps[1],yerr=errs[1], label='ZwCL2341', color='red')
    plt.errorbar(bins[1],comp_corrs[1],yerr=errs[1],linestyle='--', color='red')
    plt.xlabel(r'Log Flux $S_T$ [mJy]',size=20)
    plt.ylabel('Detected Fraction',size=20)
    plt.tick_params(axis='both', which='major', labelsize=18, length=5, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=18, length=5, width=1)
    #plt.axvline(x=sigma, color='red', linestyle='--', label=r'detection    threshold $4 \sigma$')
    plt.axhline(y=1,color='green')
    plt.xscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/ratio_plot.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
    plt.close()

if "__main__":
    names=['A2631','Zwcl2341']  
    bins=[]
    comps=[]
    comp_corrs=[]
    errs=[]
    for name in names:
        file1=name+'/'+'visib_correction_cut.txt'
        bin,comp,comp_corr,err=get_data(file1)
        bins.append(bin)
        comps.append(comp)
        errs.append(err)
        comp_corrs.append(comp_corr)
    plot(bins,comps,comp_corrs,errs)