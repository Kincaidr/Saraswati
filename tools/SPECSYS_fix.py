from astropy.io import fits

file = 'COSMOS/vla_3ghz_msmf.fits'  # Update this with your actual non-working file path

from astropy.io import fits

# Open the non-working FITS file to get the data
hdulist_non_working = fits.open(file)
non_working_data = hdulist_non_working[0].data  # Extract data from non-working file

# Create a new header
new_header = fits.Header()
new_header['SIMPLE']   = True
new_header['BITPIX']   = -32
new_header['NAXIS']    = 4
new_header['NAXIS1']   = 30000
new_header['NAXIS2']   = 30000
new_header['NAXIS3']   = 1
new_header['NAXIS4']   = 1
new_header['EXTEND']   = True
new_header['BSCALE']   = 1.0
new_header['BZERO']    = 0.0
new_header['BMAJ']     = 2.0833e-4
new_header['BMIN']     = 2.0833e-4
new_header['BPA']      = 0.0
new_header['BTYPE']    = 'Intensity'
new_header['OBJECT']   = 'COSMOS3'
new_header['BUNIT']    = 'JY/BEAM'
new_header['EQUINOX']  = 2000.0
new_header['RADESYS']  = 'FK5'
new_header['LONPOLE']  = 180.0
new_header['LATPOLE']  = 2.20583333333
new_header['PC01_01']  = 1.0
new_header['PC02_01']  = 0.0
new_header['PC03_01']  = 0.0
new_header['PC04_01']  = 0.0
new_header['PC01_02']  = 0.0
new_header['PC02_02']  = 1.0
new_header['PC03_02']  = 0.0
new_header['PC04_02']  = 0.0
new_header['PC01_03']  = 0.0
new_header['PC02_03']  = 0.0
new_header['PC03_03']  = 1.0
new_header['PC04_03']  = 0.0
new_header['PC01_04']  = 0.0
new_header['PC02_04']  = 0.0
new_header['PC03_04']  = 0.0
new_header['PC04_04']  = 1.0

# Add WCS elements for the image
new_header['CTYPE1'] = 'RA---SIN'
new_header['CRVAL1'] = 150.119166668
new_header['CDELT1'] = -5.555555617e-05
new_header['CRPIX1'] = 15000.0
new_header['CUNIT1'] = 'deg'
new_header['CTYPE2'] = 'DEC--SIN'
new_header['CRVAL2'] = 2.20583333333
new_header['CDELT2'] = 5.555555617e-05
new_header['CRPIX2'] = 15000.0
new_header['CUNIT2'] = 'deg'
new_header['CTYPE3'] = 'FREQ'
new_header['CRVAL3'] = 2.99937887379e9
new_header['CDELT3'] = 2.026758016e9
new_header['CRPIX3'] = 1.0
new_header['CUNIT3'] = 'Hz'
new_header['CTYPE4'] = 'STOKES'
new_header['CRVAL4'] = 1.0
new_header['CDELT4'] = 1.0
new_header['CRPIX4'] = 1.0
new_header['CUNIT4'] = ' '
# Add additional fields
new_header['PV2_1']   = 0.0
new_header['PV2_2']   = 0.0
new_header['TELESCOP'] = 'EVLA'
new_header['OBSERVER'] = 'Dr. Vern'
new_header['DATE-OBS'] = '2012-11-28T00:00:00.000000'
new_header['TIMESYS']  = 'UTC'
new_header['OBSRA']    = 149.535300873
new_header['OBSDEC']   = 1.62249999922
new_header['OBSGEO-X'] = -1.601156673287e6
new_header['OBSGEO-Y'] = -5.041988986066e6
new_header['OBSGEO-Z'] = 3.554879236821e6
new_header['INSTRUME'] = ' '
new_header['ALTRVAL']  = -1.52516981159e8
new_header['ALTRPIX']  = 1.0
new_header['DATE']     = '2015-07-22T08:44:01.610000'
new_header['ORIGIN']   = 'CASA 4.2.2 (prerelease r30986)'
new_header['SPECSYS']   = 'TOPOCENT' 

out_file='COSMOS/vla_3ghz_msmf.fits'
fits.writeto(out_file, data=non_working_data, header=new_header, overwrite=True)
hdulist_non_working.close()


