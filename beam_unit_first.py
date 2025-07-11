from astropy.io import fits
import sys

image = sys.argv[1]

# Define beam parameters
bmaj = 5 / 3600.0  # degrees
bmin = 5 / 3600.0  # degrees
bpa = 0.0

with fits.open(image, mode='update') as hdul:
    hdr = hdul[0].header

    # Remove all COMMENT cards
    hdr.remove('COMMENT', remove_all=True)

    # Update header values directly
    hdr['BUNIT'] = 'Jy/beam'
    hdr['BMAJ'] = bmaj
    hdr['BMIN'] = bmin
    hdr['BPA'] = bpa
    hdr['CTYPE3'] = 'Freq'
    hdr['CRVAL3'] = 1400.0e6
    hdr['CDELT3'] = 1.0e6
    hdr['CRPIX3'] = 1.0
    hdr['CTYPE4'] = 'STOKES'
    hdr['CRVAL4'] = 1.0
    hdr['CDELT4'] = 1.0
    hdr['CRPIX4'] = 1.0

    # Save changes
    hdul.flush()