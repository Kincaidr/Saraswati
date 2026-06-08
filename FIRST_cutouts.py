import requests
import bdsf
import numpy as np
from astropy.table import vstack, Table
import glob
from astropy.io import fits
from astropy.wcs import WCS

def trim_image(fits_file):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        nonzero = np.where(data != 0)

    if len(nonzero[0]) == 0:
        print("Warning: Image contains only zeros!")
        return

    ymin, ymax = np.min(nonzero[0]), np.max(nonzero[0])
    xmin, xmax = np.min(nonzero[1]), np.max(nonzero[1])
    trimmed_data = data[ymin:ymax+1, xmin:xmax+1]
    wcs = WCS(header)
    header['CRPIX1'] -= xmin
    header['CRPIX2'] -= ymin

    fits.writeto(fits_file, trimmed_data, header, overwrite=True)
    print(f"Trimmed image saved to {fits_file}")
    return(fits_file)

def stack_catalogs(catalogs, output_name):
    catalogs=glob.glob('/home/kincaid/Desktop/Saraswati_codes/A2631/images/cutouts/FIRST_srl_*.fits')
    catalog_arr=[]
    for table in catalogs:
        t=Table.read(table)
        catalog_arr.append(t)
    stacked_catalog = vstack(catalog_arr, join_type='outer')
    stacked_catalog.write(output_name, format='fits', overwrite=True)
    print(f"Stacked catalog written to {output_name}")
    return stacked_catalog

def catalog_generation(fits_image, output_name):
    img = bdsf.process_image(fits_image, rms_box=(40,40),rms_box_bright=(20,20),adaptive_thresh=150,thresh_isl=4.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = False, shapelet_do=False)
    img.write_catalog(outfile=output_name,format='fits', catalog_type='srl',clobber=True)
    print("Real catalog written")
    return(output_name)

def get_cutout(ra, dec):
        # coord = SkyCoord(ra_deg, dec_deg, unit='deg')
        # # Format RA/Dec in sexagesimal for the cutout server
        # ra_str = coord.ra.to_string(unit=u.hour, sep=':', pad=True)
        # dec_str = coord.dec.to_string(sep=':', alwayssign=True, pad=True)
        params = {
            'RA': ra,
            'Dec': dec,
            'Equinox': 'J2000',
            'ImageSize': 20,  # in arcminutes
            'ImageType': 'FITS'}
        url = 'https://third.ucllnl.org/cgi-bin/firstcutout'
        print(f"Requesting cutout for RA={ra}, Dec={dec}...")
        response = requests.get(url, params=params)
        if response.ok:
            filename =output+ f"first_cutout_{i+1}.fits"
            with open(filename, 'wb') as f:
                f.write(response.content)
        return filename

def convovle_cutout(cutout_file):
    from astropy.convolution import convolve, Gaussian2DKernel
    from astropy.io import fits
    with fits.open(cutout_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header 
    kernel_fwhm = np.sqrt(target_beam**2 - original_beam**2)
    kernel_stddev_pix = (kernel_fwhm / pix_size) / 2.3548  # FWHM to stddev
    kernel = Gaussian2DKernel(x_stddev=kernel_stddev_pix)
    convolved_data = convolve(data, kernel)
    convolved_filename = cutout_file.replace('.fits', '_convolved.fits')
    fits.writeto(convolved_filename, convolved_data,header=header, overwrite=True)
    print(f"Convolved cutout saved to {convolved_filename}")
    return(convolved_filename)
    
if __name__ == "__main__":
    target_beam=10
    original_beam=5
    pix_size=1.5
    name='Zwcl2341'
    output_stacked='/home/kincaid/Desktop/Saraswati_codes/A2631/catalogs/A2631_FIRST_stacked_srl.fits'
    output='/home/kincaid/Desktop/Saraswati_codes/A2631/images/cutouts/'
    FIRST_region_file='/home/kincaid/Desktop/Saraswati_codes/grid_coordinates.txt'
    with open(FIRST_region_file, 'r') as f:
        lines = f.readlines()
    catalogs=[]
    for i,line in enumerate(lines):
        ra,dec=line.strip().split()
        file_name=get_cutout(ra, dec)
        convolved_filename=convovle_cutout(file_name)
        convolved_filename=trim_image(convolved_filename)
        output_cat=output+f'{name}_FIRST_srl_{i}.fits'
        try:
            output_cat=catalog_generation(convolved_filename, output_cat)
        except RuntimeError:
            print(f"Error processing {convolved_filename}, skipping...")
            continue
        catalogs.append(output_cat)
    stack_catalogs(catalogs,output_stacked)
   
