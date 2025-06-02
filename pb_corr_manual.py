    
from astropy.io import fits
fits_image='Target_Abell_2631_EVLA_L_inf_1.image.tt0.fits'

data_hdu = fits.open(fits_image)[0]
data_image = data_hdu.data
header = data_hdu.header


fits_image_pb='Target_Abell_2631_EVLA_L_final.pb.tt0.fits'
data_hdu = fits.open(fits_image_pb)[0]
data_image_pb = data_hdu.data

new_data=data_image /data_image_pb
name='Target_Abell_2631_EVLA_L_inf_1.image.tt0_pbcorr.fits'
fits.writeto(name,data=new_data,header=header,overwrite=True)