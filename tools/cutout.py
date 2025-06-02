from astropy.wcs import WCS
from astropy.io import fits
import sys
from astropy.coordinates import StokesCoord


half_size = 3000  # Half the size of the cutout

# Open the FITS file
fits_image=sys.argv[1]
hdu = fits.open(fits_image)[0]
data = hdu.data
data=data[0,0,:,:]
header = hdu.header

# Read WCS and header information


w = WCS(header)
# Get the BMAJ, CRVAL1, and CRVAL2 values (Beam major axis, reference pixel coordinates)
BMAJ = header['BMAJ']
ref_x = header['CRVAL1']
ref_y = header['CRVAL2']

naxis = w.naxis

if naxis > 2:

    extra_coords = [1] * (naxis - 2)  # For simplicity, assuming 1 for any additional axes
    world_coords = [ref_x, ref_y] + extra_coords
else:
    world_coords = [ref_x, ref_y]

# Convert world coordinates (CRVAL1, CRVAL2, ...) to pixel coordinates using WCS
center_pix = w.wcs_world2pix([world_coords], 1)[0] 
center_x, center_y = int(center_pix[0]), int(center_pix[1])


# Define cutout size (center and half-size)
cutout_data = data[center_y - half_size:center_y + half_size, center_x - half_size:center_x + half_size]

# Update CRPIX1 and CRPIX2 to the new center of the cutout
header['CRPIX1'] -= (center_x - half_size)
header['CRPIX2'] -= (center_y - half_size)

# Save the cutout to a new FITS file
new_fits_image = fits_image.replace('.fits', '_cutout.fits')
fits.writeto(new_fits_image, data=cutout_data, header=header, overwrite=True)

