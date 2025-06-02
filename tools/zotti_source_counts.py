import numpy as np
import sys
import matplotlib.pyplot as plt

model_AGN='/home/kincaid/Desktop/Saraswati_codes/1d4GHz_m.txt'
AGN_data=np.loadtxt(model_AGN)
model_AGN_flux=AGN_data[:,0]
model_AGN_counts=AGN_data[:,1]

model_SFG='/home/kincaid/Desktop/Saraswati_codes/Starburst_spirals_counts_1.4GHz.txt'
SFG_data=np.loadtxt(model_SFG)
model_SFG_flux=SFG_data[:,0]
model_SFG_starburst_counts=SFG_data[:,1]
model_SFG_spiral_counts=SFG_data[:,2]

model_SFG_counts=model_SFG_starburst_counts+model_SFG_spiral_counts

#model_total_counts=model_SFG_counts+model_AGN_counts

S=np.logspace(-6,1,100)
plt.plot((10**model_AGN_flux), model_AGN_counts, label='AGN',color='red')
plt.plot((10**model_SFG_flux), 10**model_SFG_starburst_counts, label='Starburst',color='blue')
plt.plot((10**model_SFG_flux), 10**model_SFG_spiral_counts, label='Spirals',color='orange')
#plt.plot((10**model_SFG_flux), 10**model_total_counts, label='AGN+Starburst+Spirals',color='orange')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()