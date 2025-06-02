import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

filename = "output_table.txt" 
real_cat='A2631/A2631_srl.fits' 
data = np.loadtxt(filename)

ratio_old=data[:,0]
ratio=data[:,1]
flux=data[:,3]
flux_interp=data[:,4]

breakpoint()
f=interp1d(flux, ratio, bounds_error=False, fill_value="extrapolate")
y_new=f(flux_interp)
#y_new=np.append(y_new,np.ones(18))


plt.plot(flux,ratio , label='Original Data', color='blue')
plt.plot(flux_interp, y_new, label='Extrapolated Data', linestyle='dashed', color='red')
plt.xlabel('Flux mJy')
plt.ylabel('Detected fraction')
plt.legend()
plt.show()