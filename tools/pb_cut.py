from astropy.io.fits import Header
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from numpy import *
import sys
import numpy as np


def pbcor(v,rho):
 theta_b = 0.0167261 * (1.5e9/v)
 rhor = 1.18896*(rho*3.14/(180 * 60.0))/theta_b
 gain = (cos(pi*rhor)/(1.-4.*(rhor**2)))**2
 return gain

def seperation(a1,d1,a2,d2):
 c1=SkyCoord(a1*u.degree,d1*u.degree)
 c2=SkyCoord(a2*u.degree,d2*u.degree)
 sep = c1.separation(c2)
 return sep.arcmin

def pb_corr(file):
    hdu=fits.open(file)
    header = hdu[0].header
    data=hdu[0].data
    shape = (header["NAXIS1"], header["NAXIS2"])
    xx=shape[0]
    yy=shape[1]
    freq=header['CRVAL3']
    #cellsize=header['CDELT2']*60.0

    #XX, YY = np.meshgrid(np.arange(xx), np.arange(yy))
    xc=xx/2.
    yc=yy/2.
    # rho=cellsize*sqrt( (XX-xc)**2 + (YY-yc)**2)
    # seg=pbcor(freq,rho)
    # pb_corr = data/seg
    # fits.writeto(file[:-4]+'pb_corr.fits',pb_corr,header,overwrite=True)

    NAXIS1  = header['NAXIS1'] 
    NAXIS2  = header['NAXIS2']
    CDELT1  = header['CDELT1'] 
    CDELT2  = header['CDELT2'] 
    CRPIX1  = header['CRPIX1'] 
    CRPIX2  = header['CRPIX2']
    CRVAL1  = header['CRVAL1']
    CRVAL2  = header['CRVAL2']

    Hdr = fits.Header.fromstring("""\ 
    #NAXIS   =   2 
    NAXIS1  =   """+str(NAXIS1)+""" 
    NAXIS2  =   """+str(NAXIS2)+"""
    CDELT1  =   """+str(CDELT1)+"""
    CDELT2  =   """+str(CDELT2)+"""
    CRPIX1  =   """+str(CRPIX1)+"""
    CRPIX2  =   """+str(CRPIX2)+"""
    CRVAL1  =   """+str(CRVAL1)+"""
    CRVAL2  =   """+str(CRVAL2)+"""
    """, sep='\n')                  
    w = WCS(Hdr)
    yi,xi = indices((xx,yy))
    M,N=w.all_pix2world(xi, yi, 1) # convert all pix to wcs
    XC,YC=w.all_pix2world(xc, yc, 1) # convert central pix to wcs
    rho=seperation(XC,YC,M,N) # find the distance of every pix from the centre
    seg=pbcor(freq,rho)
    new_data = data/seg
    pb_corrected=file[:-4]+'pb_corr2.fits'
    fits.writeto(pb_corrected,new_data,header,overwrite=True)
    return(pb_corrected)

def cutout(fits_image):
    # Open the FITS file and get data and header
    data_hdu = fits.open(fits_image)[0]
    data_data = data_hdu.data
    data_header = data_hdu.header
    
    BMAJ = data_header['BMAJ']
    x1 = data_header['CRVAL1']
    x2 = data_header['CRVAL2']
    

    image_data = data_data[0, 0, :, :]  # Assuming the image is 4D (e.g. (Stokes, freq, y, x))
    
    center = 3000  # The center of the image
    half_size = 1500  # Half the cutout size (3000 // 2)
    
    # Perform the cutout
    cutout_data = image_data[center - half_size:center + half_size, center - half_size:center + half_size]
    
    # Update the WCS header to reflect the new cutout
    w = WCS(fits_image)
    w_cutout = w.slice((slice(center - half_size, center + half_size), 
                        slice(center - half_size, center + half_size)))
    
    # Adjust CRPIX values for the cutout
    data_header['CRPIX1'] -= (center - half_size)
    data_header['CRPIX2'] -= (center - half_size)
    
    # Save the cutout to a new FITS file
    new_fits_image = fits_image.replace('.fits', '_cutout.fits')
    fits.writeto(new_fits_image, data=cutout_data, header=data_header, overwrite=True)
    
    return cutout_data, data_header, w_cutout, BMAJ, x1, x2

if __name__ == "__main__":

    file=sys.argv[1]
    pb_corrected=pb_corr(file)
    cutout_data, cutout_header, cutout_wcs, BMAJ, x1, x2 = cutout(pb_corrected)