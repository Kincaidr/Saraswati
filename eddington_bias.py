import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.table import Table
import json
import config_file

def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def gamma_from_SCs(S, a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):

    a = np.array([a0, a1, a2, a3, a4, a5, a6])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma

def gamma_from_SCs_Mandal(S, a0=1.655, a1=-0.1150, a2=0.2272, a3=0.51788, a4=-0.449661, a5=0.160265, a6=-0.028541, a7=0.002041):

    a = np.array([a0, a1, a2, a3, a4, a5, a6, a7])
    logS = np.log10(S)
    gamma = 2.5 - np.sum([i * a[i] * logS**(i-1) for i in range(1, len(a))], axis=0)
    return gamma


def compute_eddington_bias(rec_catalog,inj_catalog, sim_plots_path, sigma):

    rec_cat=Table.read(rec_catalog)
    inj_cat=Table.read(inj_catalog)
    flux_rec = np.array(rec_cat['Total_flux']) *1e3
    flux_inj= np.array(inj_cat['Total_flux_inj']) *1e3
    flux=np.logspace(-1.5,3, 10000)
    gamma_theory=gamma_from_SCs(flux)
    plt.figure(figsize=(8, 6))
    plt.plot(flux, gamma_theory, label="Gamma vs Flux", color="blue")
    plt.xlabel("Flux Density (mJy)", fontsize=12)
    plt.ylabel("Gamma", fontsize=12)
    plt.xscale('log')
    plt.title("Gamma as a Function of Flux Density", fontsize=14)
    plt.legend()
    plt.savefig(sim_plots_path+'Eddington_gamma_plot.png')

    plt.figure(figsize=(8, 6))
    plt.scatter(flux_inj, flux_rec/flux_inj, label="Gamma vs Flux", color="blue", s=5, alpha=0.3)
    plt.xscale('log')
    plt.xlabel("Flux Injected (mJy)", fontsize=12)
    plt.ylabel("Flux Recovered/ Flux Injected", fontsize=12)
    plt.title("Flux inj vs Flux rec/Flux inj", fontsize=14)
    plt.axvline(x=5*sigma, color='red', linestyle='--', label=r'detection threshold $5 \sigma$')
    plt.legend()
    plt.savefig(sim_plots_path+'Eddington_uncorrected.png')
    
    SN=(flux_rec/10e-3)
    gamma=gamma_from_SCs_Mandal(flux_rec)
    S_true=flux_rec/2*(1+np.sqrt(1-((4*gamma)/(SN)**2)))
    plt.figure(figsize=(10, 6))
    plt.hist(np.log10(flux_inj), bins=50, alpha=0.5, label='Injected Flux ', color='orange', edgecolor='black', density=False)
    plt.hist(np.log10(flux_rec), bins=50, alpha=0.5, label='Recovered Flux ', color='blue', edgecolor='black', density=False)
    plt.hist(np.log10(S_true), bins=50, alpha=0.5, label='Recovered Flux Bias corrected', color='red', edgecolor='black', density=False)
    plt.xlabel("Log (Flux Density (mJy))", fontsize=12)
    plt.ylabel("Normalized Frequency", fontsize=12)
    plt.title("Histogram of Recovered Flux vs True Flux (Eddington Bias)", fontsize=14)
    plt.legend()
    plt.savefig(sim_plots_path+'Recovered_vs_true_flux_Eddington_Bias.png')
    rec_cat['Total_flux_Edding'] = S_true*1e-3

    mask = np.isnan(rec_cat['Total_flux_Edding']) 
    rec_cat_cleaned = rec_cat[~mask]
    rec_cat_cleaned.write(rec_catalog, format='fits', overwrite=True)

    # inj_cat_cleaned = inj_cat[~mask]
    # inj_cat_cleaned.write(inj_catalog, format='fits', overwrite=True)

    print('Catalog with Eddington bias flux column saved')

    plt.figure(figsize=(8, 6))
    plt.scatter(flux_inj, S_true/flux_inj, label="Gamma vs Flux", color="blue", s=5, alpha=0.3)
    plt.xscale('log')
    plt.axvline(x=5*sigma, color='red', linestyle='--', label=r'detection threshold $5 \sigma$')
    plt.xlabel("Flux Injected (mJy)", fontsize=12)
    plt.ylabel("Flux Recovered/ Flux Injected", fontsize=12)
    plt.title("Flux inj vs Flux rec/Flux inj", fontsize=14)
    plt.legend()
    plt.savefig(sim_plots_path+'Eddington_corrected.png')

if __name__ == "__main__":

    input=config_file.json
    config = read_config(input)
    folder=config['path']
    sim_path=config['sim_path']
    catalogs=folder+config['catalogs']
    sim_catalogs_path = sim_path+config['sim_catalogs_path']
    sim_plots_path = sim_path+config['sim_plots_path']
    merged_rec= sim_catalogs_path+config['merged_rec']
    merged_inj= sim_catalogs_path+config['merged_inj']
    sigma = config['sigma']
    compute_eddington_bias(merged_rec, merged_inj, sim_plots_path, sigma)