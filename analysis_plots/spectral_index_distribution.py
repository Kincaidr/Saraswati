import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
# example run: ipython3 spectral_index_distribution.py A2631/spectral/image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_998_srl.fits 998 A2631/spectral/image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1283_srl.fits 1283 A2631/spectral/image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1569_srl.fits 1569
# Function to read FITS catalog
def read_fits_catalog(file_path):
    """Reads a PyBDSF-generated FITS catalog and returns an Astropy Table."""
    with fits.open(file_path) as hdul:
        data = Table(hdul[1].data)
        data=data[data['S_Code']=='S']  # Assuming catalog is in the first extension
        mask=(data['Total_flux']/data['Peak_flux'] > 0.85) & (data['Total_flux']/data['Peak_flux'] < 1.15)
        data=data[mask]
        
    return data


def crossmatch_catalogs(cat1, cat2, radius=1.0):
    """
    Crossmatches two catalogs using sky coordinates.
    :param cat1: First catalog (Astropy Table).
    :param cat2: Second catalog (Astropy Table).
    :param radius: Matching radius in arcseconds.
    :return: Matched subsets of both catalogs.
    """
    coords1 = SkyCoord(cat1['RA'] * u.deg, cat1['DEC'] * u.deg)
    coords2 = SkyCoord(cat2['RA'] * u.deg, cat2['DEC'] * u.deg)

    idx, d2d, _ = coords1.match_to_catalog_sky(coords2)
    matched = d2d.arcsec < radius  # Apply matching radius condition

    return cat1[matched], cat2[idx[matched]]


def compute_spectral_index(S1, S2, nu1, nu2):
    """Computes the spectral index α using two flux densities at different frequencies."""
    valid = (S1 > 0) & (S2 > 0)  # Avoid log(0) issues
    return np.log(S1[valid] / S2[valid]) / np.log(nu1 / nu2)


# Main execution: Read multiple FITS catalogs and compute spectral indices
if len(sys.argv) < 3:
    print("Usage: python compute_spectral_index.py catalog1.fits freq1 catalog2.fits freq2 ...")
    sys.exit(1)

catalog_files = sys.argv[1::2]  # Extracts catalog filenames
frequencies = list(map(float, sys.argv[2::2]))  # Extracts frequencies

if len(catalog_files) < 2:
    print("Error: At least two catalogs are required to compute spectral indices.")
    sys.exit(1)

# Read all catalogs
catalogs = [read_fits_catalog(f) for f in catalog_files]

matched_fluxes = []

for i in range(len(catalogs) - 1):
    cat1, cat2 = crossmatch_catalogs(catalogs[i], catalogs[i + 1])
    matched_fluxes.append((cat1['Total_flux'], cat2['Total_flux'], frequencies[i], frequencies[i + 1]))

spectral_indices = []
for S1, S2, nu1, nu2 in matched_fluxes:
    alpha = compute_spectral_index(S1, S2, nu1, nu2)
    spectral_indices.append(alpha)
breakpoint()
# Flatten arrays
spectral_indices = np.concatenate(spectral_indices)
mask = (spectral_indices < 3) & (spectral_indices > -3)

plt.figure(figsize=(12, 8))
path='/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'
#plt.subplot(1, 2, 1)
plt.hist(spectral_indices, bins=30, color='blue', alpha=0.7)
plt.xlabel("Spectral Index (α)",size=20)
plt.ylabel("Number of Sources",size=20)
#plt.title("Spectral Index Distribution")
plt.axvline(np.median(spectral_indices), color='red', linestyle='dashed', label=f'Median: {np.median(spectral_indices):.2f}')
plt.legend(fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15, length=5, width=1)  # Increase size of major tick labels
plt.tick_params(axis='both', which='minor', labelsize=15, length=5, width=1)
plt.savefig(path+'Spectral_Index_Distribution.png')
plt.show()