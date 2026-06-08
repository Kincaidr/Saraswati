import numpy as np
import matplotlib.pyplot as plt

# Parameters
nu = 1.5  # GHz
theta_b = 57.5 * (nu / 1.5)**(-1)  # arcmin
theta_b_deg = theta_b / 60.0       # convert arcmin → degrees

# Angular range (degrees)
rho = np.linspace(0.001, 2.5, 2000)

# Cosine-taper beam model
x = 1.189 * rho / theta_b_deg
ab = (np.cos(np.pi * x) / (1 - 4 * x**2))**2

# Convert to dB
ab_db = 10 * np.log10(ab)

outpath='/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/'
plt.figure(figsize=(7,5))
plt.axvline(x=0.70, color='black', linestyle='--', linewidth=2, label=r'$\rho = 0.70$')                              
plt.plot(rho, ab_db, linewidth=2,label='Cosine-tapered field')
plt.xlabel(r'$\rho$ (deg)',size=14)
plt.ylabel(r'$a_b$ (dB)',size=14)
plt.ylim(-50, 1)
plt.xlim(0, 2.5)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=13)
plt.savefig(outpath+'cosine_taper_beam_pattern.png', dpi=300, bbox_inches='tight')
plt.show()
