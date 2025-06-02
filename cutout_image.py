from astropy.io import fits
import numpy as np

def cutout(fits_image):
    with fits.open(fits_image) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        image_data=data[0,0,:,:]
        breakpoint()

    ny, nx = image_data.shape
    cutout_size_pix = 1800
    half_size = cutout_size_pix // 2

    cx, cy = nx // 2, ny // 2
    x1, x2 = cx - half_size, cx + half_size
    y1, y2 = cy - half_size, cy + half_size

    new_data = np.full_like(image_data, np.nan)  # Use 0 instead if preferred
    new_data[y1:y2, x1:x2] = image_data[y1:y2, x1:x2]
    fits.writeto("cutout_image.fits", new_data, header, overwrite=True)

if __name__=="__main__":
    path='/home/kincaid/Desktop/Saraswati_codes/A2631/images/'
    fits_image=path+'image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits'
    cutout(fits_image)