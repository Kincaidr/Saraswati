from astropy.io import fits


fits_image = '/home/kincaid/Desktop/Saraswati_codes/Zwcl2341/images/mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.pbcor.fits'   
hdu = fits.open(fits_image)[0]
image_data = hdu.data*-1
header = fits.getheader(fits_image)

new_fits_image='Zwcl2341_inverted.fits'
fits.writeto(new_fits_image,data=image_data,header=header)