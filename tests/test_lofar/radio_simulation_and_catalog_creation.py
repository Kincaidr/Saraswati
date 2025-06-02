#!/usr/bin/env python
# robert.kincaid@epfl.ch

#Create radio catalog from simulated radio sources with known gaussian noise

import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import pickle
from scipy.integrate import quad, trapezoid
import bdsf
from scipy import interpolate
from astropy.wcs import WCS
from astropy.coordinates import StokesCoord
from astropy.coordinates import match_coordinates_sky
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.modeling.models import Gaussian2D
import glob
import re
import json


def catalog_generation(fits_image, name):
    
    img = bdsf.process_image(fits_image, rms_box=(40,20),rms_box_bright=(15,10),adaptive_thresh=150,thresh_isl=4.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = True) #spectralindex_do = True
    
    img.export_image(outfile=path+name+'_rms_map.fits',clobber=True,img_type='rms')
    img.export_image(outfile=path+name+'_res_map.fits',clobber=True,img_type='gaus_resid')
    img.write_catalog(outfile=path+name+'_srl.fits',format='fits', catalog_type='srl',clobber=True)

def is_too_close(new_source, sources, min_dist_arcsec,pix_size):
    min_distance=min_dist_arcsec*pix_size
    for source in sources:
        distance = np.sqrt((new_source[0] - source[0])**2 + (new_source[1] - source[1])**2)
        if distance < min_distance:
            return True
    return False

def source_counts(flux,survey_area,counts_freq,data_freq, Spectral_Index,nbins,corr=None,Range_x=None):
        
        flux = flux * (counts_freq / data_freq) ** Spectral_Index
        if Range_x is None:
                pers = np.linspace(0, 100, nbins+1)
                Range_x = np.percentile(flux, pers)
        else:
                Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1)
        centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
        difflr=np.diff(Range_x)
        hist, _ = np.histogram(flux, bins=Range_x)
        if corr ==None:
            source_tot=hist
        else:
            source_tot=hist*corr
        
        counts=norm(source_tot,difflr,survey_area,centres)
        counts_err=norm(np.sqrt(source_tot),difflr,survey_area,centres)   
        return(Range_x,source_tot,counts,counts_err,centres)


def place_sources(num_sources, xx, yy, min_distance, pix_size):
    sources = []
    attempts = 0
    max_attempts = num_sources * 10  # Limit attempts to avoid infinite loop
    
    while len(sources) < num_sources and attempts < max_attempts:
        xc = np.random.randint(xx[0], xx[1])
        yc = np.random.randint(yy[0], yy[1])
        new_source = (xc, yc)
        if not is_too_close(new_source, sources, min_distance,pix_size):
            sources.append(new_source)
        attempts += 1
    if len(sources) < num_sources:
        print(f"Warning: Only placed {len(sources)} sources out of {num_sources} after {max_attempts} attempts.")
    return sources


def generate_random_position_within_circle(max_radius, center_x, center_y):
        
        while True:
            # Randomly select pixel coordinates within a bounding box around the circle
            xc = np.random.randint(center_x - max_radius, center_x + max_radius)
            yc = np.random.randint(center_y - max_radius, center_y + max_radius)
            
            # Check if the selected point is within the circular region
            if (xc - center_x)**2 + (yc - center_y)**2 <= max_radius**2:
                return xc, yc

def plot_size(cat1,cat2,output_path,outname):

    Maj_inj=cat1['Maj_conv']*3600
    Maj_rec=cat2['Maj']*3600
    Min_rec=cat2['Min']*3600
    # mask=(Maj_rec/Maj_inj > 0.99) &  (Maj_rec/Maj_inj < 1.01)
    # #mask=(Maj_rec/Maj_inj > 1.5)
    # mask=np.where(mask)[0]
    # cat1_new=cat1[mask]
    geom_mean=np.sqrt(Maj_rec*Min_rec)
    plt.scatter(geom_mean,Maj_inj,s=7,alpha=0.7,color='blue')
    plt.plot(Maj_rec,Maj_rec,color='black')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('log Major axis rec',fontsize=12)
    plt.ylabel('log Major axis inj',fontsize=12)
    plt.title(outname[:-4], fontsize=12)
    print(outname+' plotted')
    plt.savefig(output_path+outname)
    plt.close()


def plot_flux(flux_inj,flux_rec,Dist, output_path, outname, title, xlabel='Total flux', ylabel='Peak flux'):
 
    scatter=plt.scatter(flux_rec,flux_inj,c=Dist,cmap='viridis', s=10, alpha=0.5)
    plt.colorbar(scatter,label='Distance from Phase Center')
    plt.plot(flux_rec,flux_rec,color='black')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(xlabel + ' [Jy]',fontsize=12)
    plt.ylabel(ylabel+' [Jy]',fontsize=12) 
    plt.title(title, fontsize=12)
    plt.savefig(output_path+outname)
    print(outname+' plotted')
    plt.close()


def plot_hist(flux_inj,flux_rec,nbins,output_path,outname, title ,min_size,xlabel='Total flux', ylabel='Peak flux'):
    
    plt.hist(np.log10(flux_inj),color='blue',alpha=0.8,label='injected catalog',bins=nbins)
    plt.hist( np.log10(flux_rec),color='orange',alpha=0.8,label='recovered catalog',bins=nbins)
    plt.xlabel(xlabel,size=12)
    plt.ylabel(ylabel,size=12)
    #plt.axvline(min_size, label="Major=BMAJ")
    plt.tick_params()
    plt.legend()
    plt.title(title,fontsize=12)
    print(outname+' plotted')
    plt.savefig(output_path+outname)
    plt.close()    


def extract_number(filename):
    match = re.search(r'\d+', filename)
    return int(match.group()) if match else 0

def output_catalogs(sim_catalogs_path,group):
    inj_cats_real=np.sort(glob.glob(sim_catalogs_path+group))
    inj_cats_real=   sorted( inj_cats_real, key=extract_number)
    return(inj_cats_real)


def cross_match(cat1,cat2,seperation,RA_1,DEC_1,RA_2,DEC_2,scale=None):
   
        DEC_inj=cat1[DEC_1]
        DEC_rec=cat2[DEC_2]
        
        if scale ==None:
            RA_inj=cat1[RA_1]+0
            RA_rec=cat2[RA_2]+0
        else:
            RA_inj=cat1[RA_1]+0
            RA_rec=cat2[RA_2]+360
        from astropy.units import UnitTypeError
        try:
            c1=SkyCoord(ra=RA_inj, dec=DEC_inj)
        except UnitTypeError:
            c1=SkyCoord(ra=RA_inj*u.degree, dec=DEC_inj*u.degree)
        try:
            c2=SkyCoord(ra=RA_rec, dec=DEC_rec)
        except UnitTypeError:
            c2=SkyCoord(ra=RA_rec*u.degree, dec=DEC_rec*u.degree)

        idx,d2d,_=match_coordinates_sky(c2,c1) #c2.match_to_catalog_sky(c1)
        max_sep = seperation * u.arcsec
        matches_within_4_arcsec = d2d < max_sep
        matched_indices_c1 = idx[matches_within_4_arcsec]
        cat1_matched = cat1[matched_indices_c1]
        cat2_matched = cat2[matches_within_4_arcsec]
        return(cat1_matched,cat2_matched)

def norm(source_num,difflr,survey_area,centres):
    source_tot1=source_num/ (difflr)
    source_tot2=source_tot1/ ((survey_area)*(np.pi/180)**2)
    source_norm=source_tot2*centres**(2.5)
    return(source_norm)

def ratio_source_count(flux, nbins, Range_x=None, equal_bins=False, log=False):

        if Range_x is None: Range_x = 10**(np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1))
        #if Range_x is None: Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1)
        centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
        hist, _ = np.histogram(flux, bins=Range_x)
        source_tot=hist
        return(Range_x, source_tot, centres)


def source_counts_v2(flux,survey_area,counts_freq,data_freq, Spectral_Index,nbins,corr=None,Range_x=None,equal_sources=True):
        
        # lower_bins = np.logspace(np.log10(flux_min), np.log10(flux_mid), num=lower_bins+1)
        # upper_bins = np.logspace(np.log10(flux_mid), np.log10(flux_max), num=upper_bins+1)
        # Range_x = np.concatenate((lower_bins[:-1], upper_bins))
        
        # centres = (Range_x[:-1] + Range_x[1:]) / 2.0
        # hist, _ = np.histogram(flux, bins=Range_x)

        flux = flux * (counts_freq / data_freq) ** Spectral_Index
        if Range_x is None:
            if equal_sources:
                pers = np.linspace(0, 100, nbins+1)
                Range_x = np.percentile(flux, pers)
            else:
                Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1)
        centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
        difflr=np.diff(Range_x)
        hist, _ = np.histogram(flux, bins=Range_x)
        if corr ==None:
            source_tot=hist
        else:
            source_tot=hist*corr
        
        counts=norm(source_tot,difflr,survey_area,centres)
        counts_err=norm(np.sqrt(source_tot),difflr,survey_area,centres)   
        return(Range_x,source_tot,counts,counts_err,centres)

def image_properties(fits_image):
    data_hdu = fits.open(fits_image)[0]
    data_data = data_hdu.data
    data_header = data_hdu.header
    BMAJ=data_header['BMAJ']
    BMIN=data_header['BMIN']
    BPA=data_header['BPA']
    x1=data_header['CRVAL1']
    x2=data_header['CRVAL2']
    new_data = data_data[:,:]
    w = WCS(fits_image)
    return(new_data,data_header,w, BMAJ,BMIN,BPA, x1, x2)


def gaussianv2(pix_size, A, re):
    sigma = 7
    size = int(np.round(pix_size * re)*sigma)
    x = np.arange(0, size, dtype=np.float64)[:, None]
    y = np.arange(0, size, dtype=np.float64)[None, :]
    gaussian_model = Gaussian2D(amplitude=A, x_mean=size/2, y_mean=size/2, x_stddev=re, y_stddev=re)
    gaussian_source = gaussian_model(x, y)
    return gaussian_source

def convolve_gaussian(maj_fwhm, min_fwhm, pix_size):

    sigma_x = maj_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    sigma_y = min_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    return(maj_fwhm, min_fwhm, sigma_x, sigma_y)

def gaussianv1(pix_size, A, re1,re2, theta):
    sigma = 7
    size1 = int(np.round(pix_size * re1)*sigma)
    size2 = int(np.round(pix_size * re2)*sigma)
    x = np.arange(0, size1, dtype=np.float64)[:, None]
    y = np.arange(0, size2, dtype=np.float64)[None, :]
    gaussian_model = Gaussian2D(amplitude=A, x_mean=size1/2, y_mean=size2/2, x_stddev=re1, y_stddev=re2, theta=theta)
    gaussian_source = gaussian_model(x, y)
    return gaussian_source

def SCs_Mandal(S,  a0=1.655, a1=-0.1150, a2=0.2272, a3=0.51788, a4=-0.449661, a5=0.160265, a6=-0.028541, a7=0.002041):
    # counts_freq=1400
    # data_freq=325
    # Si=-0.7
    # S = S * (counts_freq / data_freq) ** Si
    a = np.array([a0, a1, a2, a3, a4, a5, a6, a7])
    ivals = np.arange(len(a))
    vals = np.zeros_like(S)
    logS = np.log10(S)
    for i in range(len(S)):
        vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals


def SCs_Risely(S,  a0=3.192, a1=-0.223, a2=-0.846, a3=-0.261, a4=-0.024):
    a=np.array([a0, a1, a2, a3, a4])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def source_size_dist(theta_med,thetas,q=0.62,m=0.3,k=2):
    x=np.exp(-np.log(2)*(thetas/theta_med)**q)
    return(x)

def get_source_sizes(S, min_size,max_size,m=0.3, k=2, numpoints=10000, size=1,):
    theta_med=k*S**m
    thetas = np.linspace(min_size, max_size, numpoints)
    size_prob = source_size_dist(theta_med,thetas)
    size_prob_norm = size_prob/size_prob.sum()
    CDF=np.cumsum(size_prob_norm)
    inverse_CDF=interpolate.interp1d(CDF,thetas,fill_value="extrapolate")
    x=np.random.uniform(0, 1, 1)
    source_size=inverse_CDF(x)
    return source_size

def distance(x1,x2,y1,y2):
    dist=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dist

def get_source_fluxes(Svals, numpoints=100000):
    dNdSS25 = 10 ** (SCs_Bondi(Svals))
    flux_prob = dNdSS25 * (Svals/1000)**(-2.5)
    flux_prob_norm=flux_prob/np.sum(flux_prob)
    CDF=np.cumsum(flux_prob_norm)
    inverse_CDF=interpolate.interp1d(CDF,Svals,fill_value="extrapolate")
    x=np.random.uniform(0,1,numpoints)
    flux_samples=inverse_CDF(x)
    return flux_samples

def tot_num_sources(Svals, survey_area):
    dNdSvals = 10 ** (SCs_Bondi(Svals)) * (Svals/1000)**(-2.5)
    integ = trapezoid(dNdSvals, Svals/1000)
    number_sources= int(np.round(integ*(survey_area*(np.pi/180)**2)))
    return(number_sources)

def distance(x1,x2,y1,y2):
    dist=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dist

def source_properties(xc, yc, flux, BMAJ,BMIN, BPA, maj_fwhm, min_fwhm, pix_size, w, padding):
        
        Total_flux= flux
        _, _, sigma_x, sigma_y=convolve_gaussian(maj_fwhm, min_fwhm, pix_size)
        area_source= min_fwhm* maj_fwhm
        area_beam=BMAJ*BMIN
        A=Total_flux*(area_beam/area_source)   
        gaussian=gaussianv1(pix_size,A,sigma_x, sigma_y, BPA)
        x1=int(xc-(gaussian.shape[0]/2))
        x2=int(xc+(gaussian.shape[0]/2))
        y1=int(yc-(gaussian.shape[1]/2))
        y2=int(yc+(gaussian.shape[1]/2))
        xc, yc=xc-padding/2, yc-padding/2   
        aux_deg = w.pixel_to_world(yc, xc, 1,StokesCoord(1))[0]
        xc_deg, yc_deg = aux_deg.ra.value, aux_deg.dec.value  
        return(x1,x2,y1,y2,xc_deg, yc_deg, gaussian,maj_fwhm, min_fwhm,A, Total_flux)


def simulation(fits_image,lofar_cat,sim_path,sim_images_path,sim_catalogs_path, no_of_simulations) :
    output_cat_path=sim_path+sim_catalogs_path
    output_image_path=sim_path+sim_images_path
    data,header,w, BMAJ,BMIN,BPA, RA_centre, DEC_centre=image_properties(fits_image) #extract fits data, header 
    padding=2500
    paddiv2=int(padding/2)
    original_image = np.zeros((header['NAXIS1'], header['NAXIS2']))
    padded_image= np.zeros((original_image.shape[0]+padding, original_image.shape[1]+padding))
    noise=1e-9
    mean,stddev=0,noise
    uniform_noise = np.random.normal(mean, stddev, original_image.shape)
    real_noise=data
    xx=(paddiv2,paddiv2+original_image.shape[0])
    yy=(paddiv2,paddiv2+original_image.shape[1])
    pix_size=header['CDELT2']*3600
    BMAJ=BMAJ*3600
    BMIN=BMIN*3600
    number_sources=len(lofar_cat)
    print('Total number of sources', number_sources)
    flux_lofar=lofar_cat['tot_1.28GHz']
    maj_conv=lofar_cat['Maj_conv']
    min_conv=lofar_cat['Min_conv']

    for j in range(no_of_simulations):
        Data_real=np.zeros(padded_image.shape)
        RA=np.zeros(number_sources)
        DEC=np.zeros(number_sources)
        RA_pix=np.zeros(number_sources)
        DEC_pix=np.zeros(number_sources)
        Peak=np.zeros(number_sources)
        Maj_conv=np.zeros(number_sources)
        Min_conv=np.zeros(number_sources)
        Flux_total=np.zeros(number_sources)
        Dist=np.zeros(number_sources)
        sources=place_sources(number_sources, xx, yy, min_distance=15, pix_size=pix_size)
      
        for i in range(number_sources): 
            xc = sources[i][0]
            yc = sources[i][1]                                                                                      
            x1,x2,y1,y2,xc_deg, yc_deg, gaussian,_, _, A, Total_flux=source_properties(xc, yc, flux_lofar[i]*10**-3, BMAJ,BMIN,BPA, maj_conv[i], min_conv[i], pix_size,w, padding)
            Data_real[x1:x2,y1:y2] +=gaussian
            RA[i]=xc_deg
            DEC[i]=yc_deg
            RA_pix[i]=xc
            DEC_pix[i]=yc
            Peak[i]=A
            Dist[i]=distance(RA[i],RA_centre,DEC[i],DEC_centre)
            Maj_conv[i]=maj_conv[i]
            Min_conv[i]=min_conv[i]
            Flux_total[i]=Total_flux
            print('Source ' + str(i) + ' written')
        
        final_image=Data_real[paddiv2:-paddiv2,paddiv2:-paddiv2]
        sim_image_test=output_image_path+'simulated_image_test_'+str(j)+'.fits'
        sim_image_real=output_image_path+'simulated_image_real_'+str(j)+'.fits'
        fits.writeto(sim_image_real,data=final_image+real_noise,header=header,overwrite=True)
        fits.writeto(sim_image_test,data=final_image+uniform_noise,header=header,overwrite=True)
        print(sim_image_test+ ' ' +'written')
        print(sim_image_real+ ' ' +'written')
        col1 = fits.Column(name='Total_flux_inj', format='D',array=Flux_total)
        col2 = fits.Column(name='Peak_flux_inj', format='D',array=Peak)
        col3 = fits.Column(name='RA_inj', format='D',array=RA)
        col4 = fits.Column(name='DEC_inj', format='D',array=DEC)
        col5 = fits.Column(name='RA_pix', format='D',array=RA_pix)
        col6 = fits.Column(name='DEC_pix', format='D',array=DEC_pix)
        col7 = fits.Column(name='Maj_conv', format='D',array=Maj_conv)
        col8 = fits.Column(name='Min_conv', format='D',array=Min_conv)
        col9 = fits.Column(name='Dist', format='D',array=Dist)
        hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        inj_cat_real=output_cat_path+'injected_cat_'+ str(j) +'.fits'
        hdu.writeto(inj_cat_real,overwrite=True)
        print(inj_cat_real+ ' ' +'written')

def recovered_catalogs(sim_path, sim_images_path,sim_catalogs_path,no_of_simulations):

        sim_images=sim_path+sim_images_path
        sim_catalogs=sim_path+sim_catalogs_path

        test_images=np.sort(glob.glob(sim_images+'simulated_image_test_*.fits'))
        real_images=np.sort(glob.glob(sim_images+'simulated_image_real_*.fits')) 
        test_images = sorted(test_images, key=extract_number)
        real_images = sorted(real_images, key=extract_number)
        for j in range(no_of_simulations):
            img = bdsf.process_image(test_images[j], rms_box=(40,40), rms_box_bright=(15,5), adaptive_thresh=150, thresh_isl=4.0, thresh_pix=5.0,
                        interactive=False, clobber=True, spectralindex_do = False, atrous_do = False)
            
            sim_image_test=sim_catalogs+'simulated_image_test_srl_'+str(j)+'.fits'
            img.write_catalog(outfile=sim_image_test, format='fits', catalog_type='srl', clobber=True)
            print('test catalog written')  

            img = bdsf.process_image( real_images[j], rms_box=(40,40), rms_box_bright=(15,5), adaptive_thresh=150, thresh_isl=4.0, thresh_pix=5.0,
                         interactive=False, clobber=True, spectralindex_do = False, atrous_do = False)
            
            sim_image_real=sim_catalogs+'simulated_image_real_srl_'+str(j)+'.fits'
            img.write_catalog(outfile=sim_image_real, format='fits', catalog_type='srl', clobber=True)
            print('real catalog written')
        return()

def false_detections(sim_images_path, sim_catalogs_path):

    rec_images=output_catalogs(sim_images_path,group='simulated_image_real_*.fits')
    for i in range(len(rec_images)):
        #output=rec_images[i].replace('simulated_image_real_'+str(i)+'.fits','simulated_image_real_'+str(i)+'_inverted.fits')
        output=sim_images_path+'inverted_'+str(i)+'.fits'
        hdu = fits.open(rec_images[i])[0]
        data = hdu.data
        header = hdu.header
        new_data = data*-1
        fits.writeto(output,data=new_data,header=header,overwrite=True)
        print(output+ ' written')

    
def false_detections_catalogs(sim_images_path, sim_catalogs_path):
    
    images=output_catalogs(sim_images_path,group='inverted_*.fits')
    path=sim_catalogs_path 

    for i in range(len(images)):
        output=path+'simulated_image_real_'+str(i)+'_inverted_srl.fits'
        img = bdsf.process_image(images[i], rms_box=(40,40),rms_box_bright=(15,15),adaptive_thresh=150,thresh_isl=4.0,thresh_pix=5.0,
                        detection_image=images[i],interactive=False,clobber=True,spectralindex_do = False,atrous_do = False) #spectralindex_do = True
        img.write_catalog(outfile=output,format='fits', catalog_type='srl',clobber=True)
        print(output+ 'written')


def merg_catalogs(path, sim_path, sim_catalogs_path, catalogs):

    catalogs_path=sim_path+sim_catalogs_path
    output_cat=path+ catalogs
    inj_cats=output_catalogs(catalogs_path,group='injected_cat_*.fits')
    rec_cats_real=output_catalogs(catalogs_path,group='simulated_image_real_srl_*.fits')
    rec_cats_test=output_catalogs(catalogs_path,group='simulated_image_test_srl_*.fits')
    inv_cats=output_catalogs(catalogs_path,group='simulated_image_real_*_inverted_srl.fits')
    inj_cats = [Table.read(file) for file in inj_cats]
    rec_cats_real = [Table.read(file) for file in rec_cats_real]
    rec_cats_test = [Table.read(file) for file in rec_cats_test]
    inv_cats = [Table.read(file) for file in inv_cats]
    inj_real_arr=[]
    rec_test_arr=[]
    inj_test_arr=[]
    rec_real_arr=[]
    rec_real=[]
    rec_test=[]
    inj=[]
    inv_cats_arr=[]
    
    for cat in inv_cats:
        inv_cats_arr.append(cat)

    for i in range(len(rec_cats_test)) :
        rec_real.append(rec_cats_real[i])
        rec_test.append(rec_cats_test[i])
        inj.append(inj_cats[i])
    
    for i in range(len(rec_cats_test)) :
        inj_cat_test,rec_cat_test=cross_match(inj_cats[i],rec_cats_test[i],seperation=1,RA_1='RA_inj',DEC_1='DEC_inj',RA_2='RA',DEC_2='DEC',scale=360)
        inj_test_arr.append(inj_cat_test)
        rec_test_arr.append(rec_cat_test)

    for i in range(len(rec_cats_real)):
        inj_cat_real,rec_cat_real=cross_match(inj_cats[i],rec_cats_real[i],seperation=1,RA_1='RA_inj',DEC_1='DEC_inj',RA_2='RA',DEC_2='DEC',scale=360)
        inj_real_arr.append(inj_cat_real)
        rec_real_arr.append(rec_cat_real)
    
    combined_inj_real = vstack(inj_real_arr)
    combined_inj_test = vstack(inj_test_arr)
    combined_rec_real = vstack(rec_real_arr)
    combined_rec_test = vstack(rec_test_arr)
    combined_rec=vstack(rec_real)
    combined_inj=vstack(inj)
    #combined_inv=vstack(inv_cats_arr)

    combined_inj_test.write( catalogs_path+'merged_inj_test.fits', overwrite=True)
    print('merged_inj_test.fits written' )
    combined_inj_real.write( catalogs_path+'merged_inj_real.fits', overwrite=True)
    print('merged_inj_real.fits written' )
    combined_rec_test.write( catalogs_path+'merged_rec_test.fits', overwrite=True)
    print('merged_rec_test.fits written' )
    combined_rec_real.write( catalogs_path+'merged_rec_real.fits', overwrite=True)
    print('merged_rec_real.fits written' )

    combined_inj.write(catalogs_path+'merged_inj.fits', overwrite=True)
    combined_inj.write(output_cat+'merged_inj.fits', overwrite=True)
    print('merged_inj.fits written' )
    combined_rec.write(catalogs_path+'merged_rec.fits', overwrite=True)
    combined_rec.write(output_cat+'merged_rec.fits', overwrite=True)
    print('merged_rec.fits written' )
    # combined_inv.write(output_cat+'merged_inverted.fits', overwrite=True)
    # print('merged_inverted.fits' )
   

def plots(merged_inj,merged_rec,merged_inj_test,merged_inj_real,merged_rec_real,merged_rec_test,sim_path,sim_plots_path,min_size):
    
        inj_cat_test = Table.read(merged_inj_test)
        inj_cat_real = Table.read(merged_inj_real)
        rec_cat_test = Table.read(merged_rec_test)
        rec_cat_real = Table.read(merged_rec_real)
        inj_cat=Table.read(merged_inj)
        rec_cat=Table.read(merged_rec)

        flux_inj=inj_cat['Total_flux_inj']
        flux_rec=rec_cat['Total_flux']
        bmaj_inj=inj_cat_real['Maj_conv']*3600
        bmaj_rec=rec_cat_real['Maj'] *3600
        output_path=sim_path+sim_plots_path
        source_count_plot(flux_inj, counts_freq=1.4, data_freq=1.28, Spectral_Index=-0.7, output_path=output_path,name='Source_counts_inj')
        source_count_plot(flux_rec, counts_freq=1.4, data_freq=1.28, Spectral_Index=-0.7, output_path=output_path,name='Source_counts_rec')
        plot_hist(flux_inj ,flux_rec,nbins=20, output_path=output_path,outname='flux_histogram.pdf', xlabel='flux inj', title='flux on real-noise map', min_size=0, ylabel='flux rec' )
        plot_hist(bmaj_inj ,bmaj_rec,nbins=20, output_path=output_path,outname='bmaj_histogram.pdf', title='Major source size injected on real-noise map', min_size=min_size,xlabel='Major axis', ylabel='frequency')

        flux_inj=inj_cat_real['Total_flux_inj']
        flux_rec=rec_cat_real['Total_flux']
        peak_inj=inj_cat_real['Peak_flux_inj']
        peak_rec=rec_cat_real['Peak_flux']
        bmaj_inj=inj_cat_real['Maj_conv']
        bmaj_rec=rec_cat_real['Maj'] *3600
        Dist=inj_cat_real['Dist']
        
        #plot_size(inj_cat,rec_cat_real,output_path=output_path, outname='Source_size_real.pdf')
        plot_flux(flux_inj,flux_rec,Dist,output_path=output_path,outname='Total_flux_inj_vs_rec_real.pdf',title='Total flux inj vs Total flux rec', xlabel='flux inj', ylabel='flux rec')
        plot_flux(peak_inj,peak_rec,Dist,output_path=output_path,outname='Peak_flux_inj_vs_rec_real.pdf', title='Peak flux inj vs Peak flux rec',xlabel='Peak inj', ylabel='Peak rec' )
        plot_flux(bmaj_inj ,bmaj_rec,Dist,output_path=output_path,outname='bmaj_inj_vs_rec_real.pdf', title='bmaj inj vs bmaj rec',xlabel='bmaj inj', ylabel='bmaj rec' )
        
        flux_inj=inj_cat_test['Total_flux_inj']
        flux_rec=rec_cat_test['Total_flux']
        peak_inj=inj_cat_test['Peak_flux_inj']
        peak_rec=rec_cat_test['Peak_flux']
        bmaj_inj=inj_cat_test['Maj_conv'] 
        bmaj_rec=rec_cat_test['Maj'] *3600
        Dist=inj_cat_test['Dist']
        #plot_size(inj_cat,rec_cat_test,output_path=output_path,outname='Source_size_test.pdf') title='Major source size injected
        
        plot_flux(flux_inj,flux_rec,Dist,output_path=output_path,outname='Total_flux_inj_vs_rec_test.pdf',title='Total flux inj vs Total flux rec', xlabel='flux inj', ylabel='flux rec')
        plot_flux(peak_inj,peak_rec,Dist,output_path=output_path,outname='Peak_flux_inj_vs_rec_test.pdf', title='Peak flux inj vs Peak flux rec',xlabel='Peak inj', ylabel='Peak rec' )
        plot_flux(bmaj_inj ,bmaj_rec,Dist,output_path=output_path,outname='bmaj_inj_vs_rec_test.pdf', title='bmaj inj vs bmaj rec',xlabel='bmaj inj', ylabel='bmaj rec' )


def source_count_plot(flux, counts_freq, data_freq, Spectral_Index, output_path,name):
        
        _,num,counts,counts_err,bin_centre=source_counts(flux,survey_area,counts_freq,data_freq, Spectral_Index,nbins,corr=None,Range_x=None)
        #plt.errorbar(bin_centre*10**3, counts,yerr=counts_err,color='red',label=f'Simulated source counts' ,fmt='o')
        plt.scatter(bin_centre*10**3, counts,color='red',label=f'Simulated source counts')
        S=np.logspace(-1.5,3,10001)
        plt.plot(S,10 ** (SCs_Bondi(S))*(S/1000),label='Bondi 2008 1.4 GHz')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(10**-3,10**3)
        plt.xlabel('Flux S [mJy]',size=12)
        plt.ylabel('S^(2.5) * dN/dS',size=12)
        plt.title(name)
        plt.legend()
        plt.savefig(output_path+name)
        print(name+' plotted')

def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

if __name__ == "__main__":

    config_path = 'A2631.json'
    config = read_config(config_path)
    name=config['name']
    path=config['path']
    sim_path=config['sim_path']
    sigma = config['sigma']
    min_flux = config['min_flux']
    max_flux = config['max_flux']
    no_of_simulations = config['no_of_simulations']
    nbins = config['nbins']
    area=config["area"]
    catalogs=config['catalogs']
    pb_corr_image = path+config['pb_corr_image']
    residual_image = path+config['residual_image']
    sim_images_path = config['sim_images_path']
    sim_catalogs_path = config['sim_catalogs_path']
    sim_plots_path = config['sim_plots_path']
    real_catalog =  path+config['real_catalog']
    flux_corr_catalog=config['flux_corr_cat']
    incompleteness_cat= config["incompleteness_cat"]
    merged_rec_real= sim_path+sim_catalogs_path+config['merged_rec_real']
    merged_rec_test= sim_path+sim_catalogs_path+config['merged_rec_test']
    merged_inj_real= sim_path+sim_catalogs_path+config['merged_inj_real']
    merged_inj_test= sim_path+sim_catalogs_path+config['merged_inj_test']
    merged_inj= catalogs+config['merged_inj']
    merged_rec= catalogs+config['merged_rec']
    merged_rec_corr= catalogs+config['merged_rec_corr']
    _,header,w, BMAJ,BMIN,BPA,_,_=image_properties(residual_image)
    min_size=BMIN*3600 #arcsec
    max_size=3*BMIN*3600 
    survey_area=area**2
    Svals = np.linspace(min_flux, max_flux, 10000)
    
    #number_sources= tot_num_sources(Svals, survey_area)
    lofar_cat=Table.read('simul_inner.fits')
    #catalog_generation(pb_corr_image, name)
    #simulation(residual_image,lofar_cat,sim_path,sim_images_path,sim_catalogs_path, no_of_simulations) # Perform simulation and injected catalogs
    #recovered_catalogs(sim_path, sim_images_path,sim_catalogs_path,no_of_simulations) # create recovered catalogs
    #false_detections(sim_images_path, sim_catalogs_path)
    #false_detections_catalogs(sim_images_path, sim_catalogs_path)
    merg_catalogs(path, sim_path, sim_catalogs_path, catalogs)
    plots(merged_inj,merged_rec,merged_inj_test,merged_inj_real,merged_rec_real,merged_rec_test,sim_path,sim_plots_path,min_size) #create total flux and source size plots for injected vs recovered





