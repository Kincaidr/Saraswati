from astropy.io import fits


name='A2631'
fits_image = f'/home/kincaid/Desktop/Saraswati_codes/{name}/images/image_DI_Clustered.DeeperDeconv.AP4.int.restored.pbcor.fits'   
hdu = fits.open(fits_image)[0]
image_data = hdu.data*-1
header = fits.getheader(fits_image)

new_fits_image=f'{name}/images/{name}_inverted.fits'
fits.writeto(new_fits_image,data=image_data,header=header,overwrite=True)