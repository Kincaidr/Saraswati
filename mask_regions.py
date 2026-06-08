import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import matplotlib.pyplot as plt

name='Zwcl2341'  # Change this to the desired cluster name
path=f"/home/kincaid/Desktop/Saraswati_codes/{name}/flux_scale/"
image_file = path+"Zwcl2341_FIRST_large.fits.img.conv.fits"
output=path+'Zwcl2341_FIRST_large.fits.img.conv_mask.fits'
with fits.open(image_file) as hdul:
    data = hdul[0].data.squeeze()
    header = hdul[0].header
    wcs = WCS(header)
    wcs = wcs.celestial
region_file = "Zwcl_box_regions.reg"
regions = Regions.read(region_file, format='ds9')
mask = np.zeros(data.shape, dtype=bool)

for region in regions:
    try:
        pixel_region = region.to_pixel(wcs)
        region_mask = pixel_region.to_mask(mode='center')

        if region_mask is not None:
            mask_part = region_mask.to_image(data.shape)
            if mask_part is not None:
                mask |= mask_part.astype(bool)
    except Exception as e:
        print(f"Skipping region due to error: {e}")

# === Step 4: Apply the mask ===
masked_data = np.array(data)
masked_data[mask] = np.nan

# === Step 5: Save masked image ===
hdu = fits.PrimaryHDU(data=masked_data, header=header)
hdu.writeto(output, overwrite=True)
print(f"✅ Masked image saved as: {output}")
