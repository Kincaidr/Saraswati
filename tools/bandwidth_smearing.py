import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.table import Table


catalog='/home/kincaid/Desktop/Saraswati_codes/simulation/sim_catalogs/merged_rec_real.fits'
#catalog='/home/kincaid/Desktop/Saraswati_codes/A2631_srl.fits'
t=Table.read(catalog)
ra=t['RA']
dec=t['DEC']
position=SkyCoord(ra, dec, frame='icrs',unit=(u.deg,u.deg))
cluster_centre=SkyCoord(str(354.4195 ), str(0.27909), frame='icrs',unit=(u.deg,u.deg))
distance= np.abs(cluster_centre.separation(position))
mask=distance.deg < 2
t_new=t[mask]
distance=distance[mask]
flux_MeerKAT=t_new['Total_flux']
peak_MeerKAT=t_new['Peak_flux']
#flux_err_MeerKAT=t_new['E_Total_flux_inj']
Ratio=(flux_MeerKAT/peak_MeerKAT)
num_bins = 5
bin_edges = np.linspace(distance.min(), distance.max(), num_bins + 1)
bin_indices = np.digitize(distance, bin_edges)
bin_medians = [np.median(Ratio[bin_indices == i]) for i in range(1, num_bins + 1)]
print('average shift',np.mean(bin_medians))
plt.plot(bin_edges[:-1] + np.diff(bin_edges) / 2, bin_medians, color='red', marker='x', label='Median flux ratio')    
#plt.plot((bin_edges[:-1] + bin_edges[:1]) / 2, bin_medians, color='red', marker='x', label='Median flux ratio')    
print('Bin medians', bin_medians)
plt.scatter(distance,Ratio, alpha=0.6,s=10,label='Compact sources')
plt.ylim(0,5)
plt.axhline(y=1, color='black',linewidth=2)
plt.xlabel('Distance from pointing centre (deg)',fontsize=15)
plt.ylabel(r'$S_{T}/S_{P}$',fontsize=15)
plt.tight_layout()
plt.legend()
plt.savefig('smearing.pdf')
plt.show() 



