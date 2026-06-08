import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

# Function to read FITS catalog
def read_fits_catalog(file_path):
    """Reads a PyBDSF-generated FITS catalog and returns an Astropy Table."""
    with fits.open(file_path) as hdul:
        data = Table(hdul[1].data)
        data=data[data['S_Code']=='S']  # Assuming catalog is in the first extension
        mask=(data['Total_flux']/data['Peak_flux'] > 0.85) & (data['Total_flux']/data['Peak_flux'] < 1.15)
        data=data[mask]
        
    return data

# Function to compute the center (RA, DEC) from multiple catalogs
def compute_catalog_center(catalogs):
    """Computes the mean RA and DEC from all sources across multiple catalogs."""
    all_ra = np.concatenate([cat['RA'] for cat in catalogs])
    all_dec = np.concatenate([cat['DEC'] for cat in catalogs])
    return np.mean(all_ra), np.mean(all_dec)

# Function to crossmatch catalogs based on positions
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


# Function to compute spectral index between two flux measurements
def compute_spectral_index(S1, S2, nu1, nu2):
    """Computes the spectral index α using two flux densities at different frequencies."""
    valid = (S1 > 0) & (S2 > 0)  # Avoid log(0) issues
    return np.log(S1[valid] / S2[valid]) / np.log(nu1 / nu2)

# Function to calculate radial distances from a given center
def compute_radial_distance(catalog, ra_center, dec_center):
    """
    Computes the angular distance of sources from a given center (RA, DEC).
    :param catalog: Source catalog (Astropy Table).
    :param ra_center: Center RA in degrees.
    :param dec_center: Center DEC in degrees.
    :return: Array of distances in arcseconds.
    """
    coords_sources = SkyCoord(catalog['RA'] * u.deg, catalog['DEC'] * u.deg)
    center = SkyCoord(ra_center * u.deg, dec_center * u.deg)
    distances = coords_sources.separation(center).arcsec
    return distances

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

# Compute center of all sources
ra_center, dec_center = compute_catalog_center(catalogs)
print(f"Computed Center: RA = {ra_center:.6f}, DEC = {dec_center:.6f}")

# Crossmatch sources across all catalogs (pairwise)
matched_fluxes = []
radial_distances = []
for i in range(len(catalogs) - 1):
    cat1, cat2 = crossmatch_catalogs(catalogs[i], catalogs[i + 1])
    matched_fluxes.append((cat1['Total_flux'], cat2['Total_flux'], frequencies[i], frequencies[i + 1]))

    # Compute radial distances for these matched sources
    radial_distances.append(compute_radial_distance(cat1, ra_center, dec_center))

# Compute spectral indices
spectral_indices = []
for S1, S2, nu1, nu2 in matched_fluxes:
    alpha = compute_spectral_index(S1, S2, nu1, nu2)
    spectral_indices.append(alpha)

# Flatten arrays
spectral_indices = np.concatenate(spectral_indices)
radial_distances = np.concatenate(radial_distances)
mask = (spectral_indices < 3) & (spectral_indices > -3)

spectral_indices =spectral_indices[mask]
radial_distances=radial_distances[mask]

# Plot spectral index distribution
#plt.figure(figsize=(12, 5))

plots='/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'
plt.figure(figsize=(12, 8))

#plt.subplot(1, 2, 1)
plt.hist(spectral_indices, bins=30, color='blue', alpha=0.7)
plt.xlabel("Spectral Index (α)",size=22)
plt.ylabel("Number of Sources",size=22)
#plt.title("Spectral Index Distribution")
plt.axvline(np.median(spectral_indices), color='red', linestyle='dashed', label=f'Median: {np.median(spectral_indices):.2f}')
plt.legend(fontsize=30)
plt.tight_layout()
plt.savefig(plots+'spectral_index_distribution.png')
plt.show()

# # Plot spectral index vs. radial distance
# plt.subplot(1, 2, 2)
# plt.scatter(radial_distances, spectral_indices, color='green', alpha=0.6)
# plt.xlabel("Radial Distance (arcsec)")
# plt.ylabel("Spectral Index (α)")
# plt.title("Spectral Index as a Function of Radius")
# plt.axhline(np.median(spectral_indices), color='red', linestyle='dashed', label=f'Median: {np.median(spectral_indices):.2f}')
# plt.legend()
# plt.grid()






