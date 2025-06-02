import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit

def func2(x,m,c): 
    lsq=(m*x) + c 
    return lsq 

hdul=fits.open('A2631_cross_match.fits')
data=Table(hdul[1].data)

flux1=data['Total_flux_1']
flux2=data['Total_flux_2']

freq1=941
freq2=1626
alpha=[]

for i in range(len(flux1)):
        S=[flux1[i],flux2[i]]
        V=[freq1,freq2]
        V=np.array((V))
        S=np.array((S)) 

        V=np.log10(V) 
        S=np.log10(S) 
        
        coeff2, var2 = curve_fit(func2,V,S) 
        index=coeff2[0] 
        
        alpha.append(index)
breakpoint()
mean=np.mean(alpha)
print('Mean is', mean)
# pl.plot(V,S,'ro') 
# pl.plot(V,yfit,'k-') 
# print (coeff2)