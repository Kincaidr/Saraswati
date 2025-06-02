from astropy.io import fits


#fits_file = 'Zwcl2341/image_DD_beam_2poly_beam_brighptsrc5_mask.int.residual.fits'
fits_file='A2631/image_DD_beam_2poly_beam_brighptsrc3_mask_nobeam.int.residual.fits'

with fits.open(fits_file, mode='update') as hdulist:
 
    header = hdulist[0].header

    # Add or update BMAJ, BMIN, BPA, and CRVAL3 keys in the header
    header['BMAJ'] =  0.00242109470021925 # Replace with the actual beam major axis value in degrees
    header['BMIN'] =  0.00200274278221229 # Replace with the actual beam minor axis value in degrees
    header['BPA'] =  152.405629838252 # Beam position angle in degrees
    header['CRVAL3'] = 1283791015.625  # Replace with the actual value for CRVAL3
    header['FREQ'] = 1283791015.625  # Replace with the actual value for CRVAL3


    hdulist.flush()

print("BMAJ, BMIN, BPA, and CRVAL3 have been successfully added to the FITS header.")