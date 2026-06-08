
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

fits_file = 'Zwcl2341/images/image_DI_Clustered.DeeperDeconv.AP.int.restored.pbcor.fits'  # Replace with your FITS file path

with fits.open(fits_file) as hdul:
    data = hdul[0].data.squeeze()
    header = hdul[0].header
    wcs = WCS(header)


valid_mask = ~np.isnan(data)

# Find bounding box of valid region
y_inds, x_inds = np.where(valid_mask)
y_min, y_max = y_inds.min(), y_inds.max()
x_min, x_max = x_inds.min(), x_inds.max()

# Trim data
trimmed_data = data[ y_min:y_max+1, x_min:x_max+1]

# Update WCS to reflect the new reference pixel location
new_wcs = wcs.deepcopy()
new_wcs.wcs.crpix[0] -= x_min  # CRPIX1
new_wcs.wcs.crpix[1] -= y_min  # CRPIX2

# Generate new header
new_header = new_wcs.to_header()
for key in ['BUNIT', 'OBJECT', 'BMAJ', 'BMIN', 'BPA']:
    if key in header:
        new_header[key] = header[key]

output_file = 'trimmed_fits.fits'
hdu = fits.PrimaryHDU(data=trimmed_data, header=new_header)
hdu.writeto(output_file, overwrite=True)
print(f"Trimmed FITS written to {output_file}")