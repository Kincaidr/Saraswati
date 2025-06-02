import numpy as np
from astropy.table import Table
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/A2631/catalogs/A2631_srl.fits')

fluxes=cat['Total_flux']*1e3

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



# Step 1: Define the hybrid log-normal + power-law model
def hybrid_model(S, A, B, mu, sigma, alpha):
    """Hybrid model combining log-normal and power-law."""
    log_normal = A * (1 / (S * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(S) - mu)**2) / (2 * sigma**2))
    power_law = B * (S**-alpha)
    return log_normal + power_law

# Step 2: Bin the data for comparison
bins = np.logspace(np.log10(min(fluxes)), np.log10(max(fluxes)), 50)  # Log-spaced bins
bin_centers = np.sqrt(bins[:-1] * bins[1:])
hist, _ = np.histogram(fluxes, bins=bins, density=True)  # Normalized histogram

# Step 3: Fit the hybrid model to the data
popt, pcov = curve_fit(hybrid_model, bin_centers, hist, p0=[1, 1, 0, 1, 2])  # Initial guesses

# Extract fitted parameters
A_fit, B_fit, mu_fit, sigma_fit, alpha_fit = popt
print(f"Fitted Parameters: A={A_fit:.2f}, B={B_fit:.2f}, mu={mu_fit:.2f}, sigma={sigma_fit:.2f}, alpha={alpha_fit:.2f}")

# Step 4: Plot the results
plt.figure(figsize=(8, 6))
plt.hist(fluxes, bins=bins, density=True, alpha=0.5, label="Data")  # Plot data histogram
S_plot = np.linspace(min(fluxes), max(fluxes), 1000)  # Smooth flux range for plotting
plt.plot(S_plot, hybrid_model(S_plot, *popt), label="Fitted Model", color="red")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Flux (mJy)")
plt.ylabel("Density")
plt.legend()
plt.show()
