


import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions, PixCoord
from astropy.coordinates import SkyCoord

# Load FITS file

fits_filename = "/home/kincaid/Desktop/Saraswati_codes/A2631/A2631_rms_map.fits"
region_filename = "/home/kincaid/Desktop/Saraswati_codes/A2631/ds9.reg"

hdul = fits.open(fits_filename)
data = hdul[0].data[0,0,:,:]  # Assuming 2D image in primary HDU
header = hdul[0].header
wcs = WCS(header)  # Load WCS for coordinate transformations

regions = Regions.read(region_filename, format="ds9")
mask = np.zeros_like(data, dtype=bool)

for region in regions:
    if hasattr(region, "center"):  # Check if it has a SkyCoord center
        skycoord = region.center   
        x, y = wcs.celestial.world_to_pixel(skycoord)
        # Convert to pixel coordinates
        pixcoord = PixCoord(x=x, y=y)
        region = region.to_pixel(wcs.celestial)
        region_mask = region.to_mask(mode="center").to_image(data.shape)
        mask |= region_mask.astype(bool)


# Apply mask to the image
data[mask] = np.nan    
fits.writeto("masked_image.fits", data, header, overwrite=True)
print("Mask applied and saved to masked_image.fits")
