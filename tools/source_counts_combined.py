import numpy as np
import sys
import matplotlib.pyplot as plt

def SCs_Mandal(S,  a0=1.655, a1=-0.1150, a2=0.2272, a3=0.51788, a4=-0.449661, a5=0.160265, a6=-0.028541, a7=0.002041):
    # counts_freq=1400
    # data_freq=325
    # Si=-0.7
    # S = S * (counts_freq / data_freq) ** Si
    a = np.array([a0, a1, a2, a3, a4, a5, a6, a7])
    ivals = np.arange(len(a))
    vals = np.zeros_like(S)
    logS = np.log10(S)
    for i in range(len(S)):
        vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals
# Load the data from the text file

def counts_plot(loaded_data):
    bin_centre = loaded_data[:, 0]
    counts = loaded_data[:, 1]
    counts_err= loaded_data[:, 2]
    counts_corr= loaded_data[:, 3]
    counts_corr_err= loaded_data[:, 4]
    return(bin_centre, counts, counts_err, counts_corr, counts_corr_err)

def flux_scale(flux,counts_freq,data_freq,si):
    flux = flux * (counts_freq / data_freq) ** si
    return(flux)

#zwcl='/home/kincaid/Desktop/Saraswati_codes/Zwcl2341/catalogs/source_counts.txt'
#zwcl_data = np.loadtxt(zwcl)

A2631='/home/kincaid/Desktop/Saraswati_codes/A2631/catalogs/source_counts.txt'
A2631_data = np.loadtxt(A2631)

model_AGN='/home/kincaid/Desktop/Saraswati_codes/1d4GHz_m.txt'
AGN_data=np.loadtxt(model_AGN)
model_AGN_flux=AGN_data[:,0]
model_AGN_counts=AGN_data[:,1]

model_SFG='/home/kincaid/Desktop/Saraswati_codes/Starburst_spirals_counts_1.4GHz.txt'
SFG_data=np.loadtxt(model_SFG)
model_SFG_flux=SFG_data[:,0]
model_SFG_starburst_counts=SFG_data[:,1]
model_SFG_spiral_counts=SFG_data[:,2]
model_SFG=model_SFG_starburst_counts+model_SFG_spiral_counts
#model_combined=model_SFG_spiral_counts+model_AGN_counts+model_SFG_spiral_counts
#zwcl_bin_centre, zwcl_counts, zwcl_counts_err,_,_=counts_plot(zwcl_data)
A2631_bin_centre, A2631_counts, A2631_counts_err, A2631_counts_corr,A2631_counts_corr_err=counts_plot(A2631_data)

# Extract the columns back into separate variables

Risely_flux=np.array([0.242,0.284,0.333,0.391,0.460,0.540,0.634,0.754,0.875,1.10,1.49,2,2.70,3.65,4.92,6.63,9.70,15.41,24.36,38.61,61.20,97,153.7,243.6,737])
Risely_counts=[12.49,12.26,14.33,15.48,15.24,15.72,15.73,18.37,18.88,22.55,22.58,28.09,31.61,41.89,52.37,58.26,82.08,168.5,215.2,251.1,325.8,520,908,1164.3,1751.2]
Resiley_error=[0.42,0.46,0.58,0.68,0.76,0.87,0.97,1.18,1.34,1.27,1.69,2.18,2.86,4.09,5.68,7.46,9.42,19,30.1,45.8,72.9,130,242.7,388.1,528]
Solmic_flux=np.array([3.67, 5.51, 8.27, 12.40, 18.60, 27.90, 41.85, 62.78, 94.17, 141.25, 211.88])
Solmic_counts=[23.04,32.24,49.17,68.25,111.83,108.73,177.52,285.36,524.18,825.42,758.19]
Solmic_error=[3.91,5.97,9.83,15.66,27.12,36.34,62.76,107.85,198.12,336.98,437.74]
owen_flux=np.array([375,475,600,900,1350,2000,3000,4500,6750,10000,20000,40000,80000])
owen_counts=[25.4,25.9,16.9,17.2,19.3,27.3,30.8,58,48.7,89,243.6,495.6,398.1]

Risely_flux=flux_scale(Risely_flux,counts_freq=1400,data_freq=325,si=-0.7)
Solmic_flux=flux_scale(Solmic_flux,counts_freq=1400,data_freq=324,si=-0.7)
MeerKAT_flux=flux_scale(A2631_bin_centre*10**3,counts_freq=1400,data_freq=1280,si=-0.7)
plt.errorbar(Risely_flux,Risely_counts, yerr=Resiley_error,color='green', label='GMRT 325 MHz Risely 2017 SuperCLASS',fmt='v')
plt.errorbar(Solmic_flux,Solmic_counts, yerr=Solmic_error,color='orange', label='GMRT 325 MHz Solmic 2014 COSMOS',fmt='^')
#plt.errorbar(zwcl_bin_centre*10**3, zwcl_counts,yerr=zwcl_counts_err,color='red',label=f'MeerKAT 1283 MHz Zwcl ' ,fmt='o')
plt.errorbar(MeerKAT_flux, A2631_counts,yerr=A2631_counts_err,color='blue',label=f'MeerKAT 1283 MHz 2019 A2631' ,fmt='+')
plt.errorbar(MeerKAT_flux, A2631_counts_corr,yerr=A2631_counts_corr_err,color='red',label=f'MeerKAT 1283 MHz 2019 A2631 corrected' ,fmt='+')
plt.errorbar(owen_flux*10**-3,owen_counts,color='cyan',label=f' VLA 325 MHz Owen 2012 SWIRE' ,fmt='d')
plt.plot((10**model_AGN_flux)*10**3, model_AGN_counts, label='AGN',color='red')
plt.plot((10**model_SFG_flux)*10**3, 10**model_SFG_starburst_counts, label='Starburst',color='blue')
plt.plot((10**model_SFG_flux)*10**3, 10**model_SFG_spiral_counts, label='Spirals',color='orange')
S=np.logspace(-1.5,3,10001)
S=flux_scale(S,counts_freq=1400,data_freq=150,si=-0.8)
#plt.plot(S,10 ** (SCs_Mandal(S)),label='Mandal 2021 150 MHz')
#plt.plot(S,10 ** (SCs_Risely(S)),label='Risely 2016 1.4 GHz')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-3,10**3)
plt.xlabel('Flux S [mJy]',size=12)
plt.ylabel('S^(2.5) * dN/dS',size=12)
plt.legend()
plt.show()
plt.close()



