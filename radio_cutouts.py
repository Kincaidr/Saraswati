import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from matplotlib.gridspec import GridSpec


def radio_cutouts(fits_image1, fits_image2):
    # cutouts = [
    #     (fits_image1,353.3307451, -0.4450689, 0.08, 0.08, 50,  1e-05, 10e-05,5,2,1),
    #     (fits_image1,354.4408546, -0.7057400, 0.05, 0.05, 30, 1e-05, 10e-05, 5,2,2),
    #     (fits_image1,354.6049207, -0.1103067, 0.03, 0.03, 30, 1e-05, 10e-05, 5,2,3),
    #     (fits_image2,356.8094313,  -0.3192631, 0.04, 0.04, 50, 1e-05, 9e-05, 5,2,4),
    #     (fits_image1,353.7103476, -0.1737802, 0.05, 0.05, 30, 1e-05, 10e-05, 5,2,5),
    #     (fits_image1,353.9070295,  0.2746370, 0.07, 0.07, 50, 1e-05, 10e-05, 5,2,6),
    #     (fits_image2,356.1199781,  0.6516063, 0.05, 0.05, 30, 1e-05, 10e-05, 5,2,7),
    #     (fits_image2,356.9196948,-0.6843053, 0.05, 0.05, 50, 1e-05, 3e-05,   5,2,8),
    #     (fits_image2, 355.8568747,0.1695829, 0.03,0.03 ,  1, 0.5e-05 ,10e-05,5,2,5)
    # ]

    vmin,vmax=0.5e-05 ,5e-05
    cutouts = [
        (fits_image1,353.3307451, -0.4450689, 0.08, 0.08, 50,5,2,1),
        (fits_image1,354.4408546, -0.7057400, 0.05, 0.05, 50, 5,2,2),
        (fits_image1,354.6049207, -0.1103067, 0.03, 0.03, 30, 5,2,3),
        (fits_image2,356.8094313,  -0.3192631, 0.04, 0.04, 50, 5,2,4),
        (fits_image2,356.1199781,  0.6516063, 0.05, 0.05, 30, 5,2,5),
        (fits_image1,353.9070295,  0.2746370, 0.07, 0.07, 50, 5,2,6),
        (fits_image2, 355.8568747,0.1695829, 0.03,0.03 ,  50,5,2,7),
        (fits_image1, 354.2668945,0.5494581, 0.025,0.025 ,  30,5,2,8),
        (fits_image1, 354.2403570, 0.8818233, 0.028,0.028 ,  30, 5,2,9),
        (fits_image1, 354.2638000, 0.9437835, 0.03,0.03 ,  30, 5,2,10)
    ]
    
    fig = plt.figure(figsize=(6, 12))     
    for i, (fits_image,ra, dec, width, height, scalebar_arcsec, r,c,p) in enumerate(cutouts):
        print(r,c,p)
        ax = plt.subplot(r,c,p)
        with fits.open(fits_image) as hdul:
            image_data = hdul[0].data
            image_data=image_data[0,0,:,:]
            wcs = WCS(hdul[0].header)

        center = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='fk5')
        x_center, y_center = skycoord_to_pixel(center, wcs)
        pixel_scale = np.abs(wcs.pixel_scale_matrix[1, 1])  # Degrees per pixel
        width_pix = width / pixel_scale
        height_pix = height / pixel_scale

        x_min, x_max = int(x_center - width_pix / 2), int(x_center + width_pix / 2)
        y_min, y_max = int(y_center - height_pix / 2), int(y_center + height_pix / 2)
        cutout = image_data[y_min:y_max, x_min:x_max]

        im = ax.imshow(cutout, cmap='gray', origin='lower', norm=plt.Normalize(vmin=vmin, vmax=vmax))
        scalebar_pix = (scalebar_arcsec / 3600) / pixel_scale  # Convert arcsec to pixels
        ax.plot([10, 8 + scalebar_pix], [12, 12], color='yellow', lw=3)
        ax.text(10 + scalebar_pix / 2, 15, f"{scalebar_arcsec}''", color='yellow', fontsize=15, ha='center')
        ax.set_xticks([])
        ax.set_yticks([])
        
    plt.subplots_adjust(wspace=0.05, hspace=0.05, right=0.85, bottom=0.05)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/radio_cutouts.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()


if "__main__":
    fits_image_A2631='/home/kincaid/Desktop/Saraswati_codes/A2631/images/mypipelinerun_ABELL2631_4-MFS-image.fits'
    fits_image_Zwcl2341='/home/kincaid/Desktop/Saraswati_codes/Zwcl2341/images/mypipelinerun_ZwCl2341_1_p_0000_4-MFS-image.fits'
    output_path='/home/kincaid/Desktop/Saraswati_codes/'
    #radio_cutout(output_path,fits_image_A2631)
    radio_cutouts(fits_image_A2631, fits_image_Zwcl2341)