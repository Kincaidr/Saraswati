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


def is_too_close(new_source, sources, min_dist_arcsec,pix_size):
    min_distance=min_dist_arcsec*pix_size
    for source in sources:
        distance = np.sqrt((new_source[0] - source[0])**2 + (new_source[1] - source[1])**2)
        if distance < min_distance:
            return True
    return False


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

def plot_flux(flux_rec, flux_inj, Dist, output, xlabel, ylabel):
        scatter=plt.scatter(flux_rec,flux_inj,c=Dist,cmap='viridis', s=10, alpha=0.5)
        plt.colorbar(scatter,label='Distance from Phase Center')
        plt.plot(flux_rec,flux_rec,color='black')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel(xlabel+' [Jy]',fontsize=12)
        plt.ylabel(ylabel+' [Jy]',fontsize=12) 
        plt.savefig(output)
        print(output+' plotted')
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

def correct_axes(fits_image):
    data_hdu = fits.open(fits_image)[0]
    data_data = data_hdu.data
    data_header = data_hdu.header
    BMAJ=data_header['BMAJ']
    data_header['BMIN']=BMAJ
    data_header['CRVAL3']=1.4e9
    new_data = data_data[0,0,:,:]
    w = WCS(fits_image)
    #new_fits_image=fits_image.replace(fits_image[-4:],'corr.fits')
    #fits.writeto(new_fits_image,data=new_data,header=data_header,overwrite=True)
    return(new_data,data_header,w, BMAJ, x1, x2)

def gaussianv2(pix_size, A, re):
    sigma = 7
    size = int(np.round(pix_size * re)*sigma)
    x = np.arange(0, size, dtype=np.float64)[:, None]
    y = np.arange(0, size, dtype=np.float64)[None, :]
    gaussian_model = Gaussian2D(amplitude=A, x_mean=size/2, y_mean=size/2, x_stddev=re, y_stddev=re)
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

def get_source_fluxes(Svals, numpoints=10000):
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

def source_properties(xc, yc, flux, BMAJ, min_size, max_size, pix_size, convert, w, padding, source_type):
        
        Total_flux= flux
        if source_type=='resolved':
            MAJ=get_source_sizes( Total_flux,min_size=min_size,max_size=max_size)
            MIN=MAJ
            re1_sigma=MAJ/convert
            re1_sigma/=pix_size
            area_beam=2*np.pi*BMAJ**2
            area_source=2*np.pi*MAJ**2
            A=Total_flux*(area_beam/area_source)   
        else:   
            MAJ=get_source_sizes( Total_flux,min_size=min_size,max_size=max_size)
            MIN=MAJ
            re1_sigma=MAJ/convert
            re1_sigma/=pix_size
            A=Total_flux

        gaussian=gaussianv2(pix_size,A,re1_sigma)
        x1=int(xc-(gaussian.shape[0]/2))
        x2=int(xc+(gaussian.shape[0]/2))
        y1=int(yc-(gaussian.shape[1]/2))
        y2=int(yc+(gaussian.shape[1]/2))
        xc, yc=xc-padding/2, yc-padding/2   
        aux_deg = w.pixel_to_world(yc, xc, 1,StokesCoord(1))[0]
        xc_deg, yc_deg = aux_deg.ra.value, aux_deg.dec.value  
        return(x1,x2,y1,y2,xc_deg, yc_deg, gaussian, MAJ, MIN,A, Total_flux)

def simulation(fits_image,pb_fits_image, min_size,max_size,sim_images_path,no_of_simulations, number_sources, Svals) :

    data,header,w, BMAJ, RA_centre, DEC_centre=correct_axes(fits_image) #extract fits data, header 
    data_pb,header_pb,w_pb, BMAJ, x2 ,y2=correct_axes(pb_fits_image)
    padding=2500
    paddiv2=int(padding/2)
    original_image = np.zeros((header['NAXIS1'], header['NAXIS2']))
    padded_image= np.zeros((original_image.shape[0]+padding, original_image.shape[1]+padding))

    pb_noise = data_pb
    real_noise=data
    xx=(paddiv2,paddiv2+original_image.shape[0])
    yy=(paddiv2,paddiv2+original_image.shape[1])
    print('Total number of sources', number_sources)
    pix_size=header['CDELT2']*3600
    BMAJ=BMAJ*3600
    convert = 2.0 * np.sqrt(2.0 * np.log(2))
    flux_samples=get_source_fluxes(Svals)

    for j in range(no_of_simulations):
        Data_real=np.zeros(padded_image.shape)
        RA=np.zeros(number_sources)
        DEC=np.zeros(number_sources)
        RA_pix=np.zeros(number_sources)
        DEC_pix=np.zeros(number_sources)
        Peak=np.zeros(number_sources)
        Maj=np.zeros(number_sources)
        Min=np.zeros(number_sources)
        Flux_total=np.zeros(number_sources)
        Dist=np.zeros(number_sources)
        sources=place_sources(number_sources, xx, yy, min_distance=15, pix_size=pix_size)
      
        for i in range(number_sources): 
            xc = sources[i][0]
            yc = sources[i][1]
            x1,x2,y1,y2,xc_deg, yc_deg, gaussian, MAJ,MIN, A, Total_flux=source_properties(xc, yc, flux_samples[i]*10**-3, BMAJ, min_size, max_size, pix_size, convert,w, padding, source_type='resolved')
            Data_real[x1:x2,y1:y2] +=gaussian
            RA[i]=xc_deg
            DEC[i]=yc_deg
            RA_pix[i]=xc
            DEC_pix[i]=yc
            Peak[i]=A
            Dist[i]=distance(RA[i],RA_centre,DEC[i],DEC_centre)
            Maj[i]=MAJ/3600
            Min[i]=MIN/3600
            Flux_total[i]=Total_flux
            print('Source ' + str(i) + ' written')

        final_image=Data_real[paddiv2:-paddiv2,paddiv2:-paddiv2]
        sim_image_test=sim_images_path+'simulated_image_test_'+str(j)+'.fits'
        sim_image_real=sim_images_path+'simulated_image_real_'+str(j)+'.fits'
        fits.writeto(sim_image_real,data=final_image+real_noise,header=header,overwrite=True)
        fits.writeto(sim_image_test,data=final_image+pb_noise,header=header,overwrite=True)
        print(sim_image_test+ ' ' +'written')
        print(sim_image_real+ ' ' +'written')
        col1 = fits.Column(name='Total_flux_inj', format='D',array=Flux_total)
        col2 = fits.Column(name='Peak_flux_inj', format='D',array=Peak)
        col3 = fits.Column(name='RA_inj', format='D',array=RA)
        col4 = fits.Column(name='DEC_inj', format='D',array=DEC)
        col5 = fits.Column(name='RA_pix', format='D',array=RA_pix)
        col6 = fits.Column(name='DEC_pix', format='D',array=DEC_pix)
        col7 = fits.Column(name='Maj_inj', format='D',array=Maj)
        col8 = fits.Column(name='Min_inj', format='D',array=Min)
        col9 = fits.Column(name='Dist', format='D',array=Dist)
        hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        inj_cat_real=sim_catalogs_path+'injected_cat_'+ str(j) +'.fits'
        hdu.writeto(inj_cat_real,overwrite=True)
        print(inj_cat_real+ ' ' +'written')

def recovered_catalogs(sim_images_path,sim_catalogs_path,no_of_simulations):

        test_images=np.sort(glob.glob(sim_images_path+'simulated_image_test_*.fits'))
        real_images=np.sort(glob.glob(sim_images_path+'simulated_image_real_*.fits')) 
        test_images = sorted(test_images, key=extract_number)
        real_images = sorted(real_images, key=extract_number)
        
        for j in range(no_of_simulations):
            
            img = bdsf.process_image(test_images[j], rms_box=(40,40), rms_box_bright=(15,15), adaptive_thresh=150, thresh_isl=4.0, thresh_pix=5.0,
                        interactive=False, clobber=True, spectralindex_do = False, atrous_do = False)
            sim_image_test=sim_catalogs_path+'simulated_image_real_pb_srl_'+str(j)+'.fits'
            img.write_catalog(outfile=sim_image_test, format='fits', catalog_type='srl', clobber=True)
            print('test catalog written')  
            img = bdsf.process_image( real_images[j], rms_box=(40,40), rms_box_bright=(15,15), adaptive_thresh=150, thresh_isl=4.0, thresh_pix=5.0,
                         interactive=False, clobber=True, spectralindex_do = False, atrous_do = False)
            sim_image_real=sim_catalogs_path+'simulated_image_real_srl_'+str(j)+'.fits'
            img.write_catalog(outfile=sim_image_real, format='fits', catalog_type='srl', clobber=True)
            print('real catalog written')
        return()

def merg_catalogs(sim_catalogs_path):
    
    inj_cats=output_catalogs(sim_catalogs_path,group='injected_cat_*.fits')
    rec_cats_real=output_catalogs(sim_catalogs_path,group='simulated_image_real_pb_srl_*.fits')
    rec_cats_test=output_catalogs(sim_catalogs_path,group='simulated_image_real_srl_*.fits')
    inj_cats = [Table.read(file) for file in inj_cats]
    rec_cats_pb = [Table.read(file) for file in rec_cats_real]
    rec_cats_nopb = [Table.read(file) for file in rec_cats_test]
    rec_cats_pb_arr=[]
    rec_cats_nopb_arr=[]
    inj_pb_arr=[]
    rec_pb_arr=[]
    inj_nopb_arr=[]
    rec_nopb_arr=[]

    for i in range(len(rec_cats_test)) :
        rec_cat_pb,rec_cat_nopb=cross_match(rec_cats_pb[i],rec_cats_nopb[i],seperation=1,RA_1='RA',DEC_1='DEC',RA_2='RA',DEC_2='DEC',scale=360)
        rec_cats_pb_arr.append(rec_cat_pb)
        rec_cats_nopb_arr.append(rec_cat_nopb)
  
    for i in range(len(rec_cats_test)) :
        inj_cat_pb, rec_cat_pb=cross_match(inj_cats[i],rec_cats_pb[i],seperation=1,RA_1='RA_inj',DEC_1='DEC_inj',RA_2='RA',DEC_2='DEC',scale=360)
        inj_pb_arr.append(inj_cat_pb)
        rec_pb_arr.append(rec_cat_pb)

    for i in range(len(rec_cats_test)) :
        inj_cat_nopb, rec_cat_nopb=cross_match(inj_cats[i],rec_cats_nopb[i],seperation=1,RA_1='RA_inj',DEC_1='DEC_inj',RA_2='RA',DEC_2='DEC',scale=360)
        inj_nopb_arr.append(inj_cat_nopb)
        rec_nopb_arr.append(rec_cat_nopb)

    combined_rec_nopb = vstack(rec_cats_nopb_arr)
    combined_rec_pb = vstack(rec_cats_pb_arr)

    inj_pb = vstack(inj_pb_arr)
    rec_pb = vstack(rec_pb_arr)
    inj_nopb = vstack(inj_nopb_arr)
    rec_nopb = vstack(rec_nopb_arr)

    breakpoint()

    inj_pb.write(sim_catalogs_path+'inj_pb.fits', overwrite=True)
    print('inj_pb.fits' )
    rec_pb.write(sim_catalogs_path+'rec_pb.fits', overwrite=True)
    print('rec_pb.fits' )

    inj_nopb.write(sim_catalogs_path+'inj_nopb.fits', overwrite=True)
    print('inj_nopb.fits' )
    rec_nopb.write(sim_catalogs_path+'rec_nopb.fits', overwrite=True)
    print('rec_nopb.fits' )

    combined_rec_nopb.write(sim_catalogs_path+'merged_rec_nopb.fits', overwrite=True)
    print('merged_rec_nopb.fits' )
    combined_rec_pb.write(sim_catalogs_path+'merged_rec_pb.fits', overwrite=True)
    print('merged_rec_pb.fits' )


def plots(merged_rec_nopb,merged_rec_pb,inj_nopb,inj_pb,rec_nopb,rec_pb,sim_plots_path):

        rec_cat_test = Table.read(merged_rec_nopb)
        rec_cat_real = Table.read(merged_rec_pb)

        inj_nopb_cat = Table.read(inj_nopb)
        inj_pb_cat = Table.read(inj_pb)

        rec_nopb_cat = Table.read(rec_nopb)
        rec_pb_cat = Table.read(rec_pb)

        RA=rec_cat_test['RA']
        DEC=rec_cat_test['DEC']

        RA_rec_pb=rec_pb_cat['RA']
        DEC_rec_pb=rec_pb_cat['DEC']

        RA_rec_nopb=rec_nopb_cat['RA']
        DEC_rec_nopb=rec_nopb_cat['DEC']

        CRVAL1=-5.58083333333331
        CRVAL2=0.276666666666667
        Dist=distance(RA,CRVAL1,DEC,CRVAL2)
        sorted_indices = np.argsort(Dist)
        rec_cat_test=rec_cat_test[sorted_indices]
        rec_cat_real=rec_cat_real[sorted_indices]
        flux_rec_test= rec_cat_test['Total_flux']
        flux_rec_real= rec_cat_real['Total_flux']

        Dist_inj_pb=distance(RA_rec_pb,CRVAL1,DEC_rec_pb,CRVAL2)
        sorted_indices = np.argsort(Dist_inj_pb)
        inj_pb_cat=inj_pb_cat[sorted_indices]
        rec_pb_cat=rec_pb_cat[sorted_indices]
        flux_inj_pb_cat=inj_pb_cat['Total_flux_inj']
        flux_rec_pb_cat=rec_pb_cat['Total_flux']

        Dist_inj_nopb=distance(RA_rec_nopb,CRVAL1,DEC_rec_nopb,CRVAL2)
        sorted_indices = np.argsort(Dist_inj_nopb)
        inj_nopb_cat=inj_nopb_cat[sorted_indices]
        rec_nopb_cat=rec_nopb_cat[sorted_indices]
        flux_inj_nopb_cat=inj_nopb_cat['Total_flux_inj']
        flux_rec_nopb_cat=rec_nopb_cat['Total_flux']

        breakpoint()
        # Dist_rec_pb=distance(RA_rec_pb,CRVAL1,DEC_rec_pb,CRVAL2)
        # Dist_rec_nopb=distance(RA_rec_nopb,CRVAL1,DEC_rec_nopb,CRVAL2)
        
        plot_flux(flux_rec_real,flux_rec_test, Dist=Dist, output=sim_plots_path+'rec_pb_vs_rec_nopb.pdf', xlabel='Flux rec pb', ylabel='Flux rec no pb')
        plot_flux(flux_inj_pb_cat,flux_rec_pb_cat, Dist=Dist_inj_pb, output=sim_plots_path+'inj_vs_rec_pb.pdf' , xlabel='Flux inj pb', ylabel='Flux rec pb')
        plot_flux(flux_inj_nopb_cat,flux_rec_nopb_cat, Dist=Dist_inj_nopb, output=sim_plots_path+'inj_vs_rec_nopb.pdf' , xlabel='Flux inj nopb', ylabel='Flux rec nopb')

          

def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

if __name__ == "__main__":

    config_path = 'A2631.json'
    config = read_config(config_path)
    name=config['name']
    path=config['path']
    catalogs=config['catalogs']
    pb_residual_image = path+config['pb_residual_image']
    residual_image = path+config['residual_image']
    sigma = config['sigma']
    min_flux = config['min_flux']
    max_flux = config['max_flux']
    no_of_simulations = config['no_of_simulations']
    nbins = config['nbins']
    area=config["area"]
    sim_images_path = path+config['sim_images_path']
    sim_catalogs_path = path+config['sim_catalogs_path']
    sim_plots_path = path+config['sim_plots_path']
    merged_rec_pb= sim_catalogs_path+config['merged_rec_pb']
    merged_rec_nopb= sim_catalogs_path+config['merged_rec_nopb']
    inj_nopb= sim_catalogs_path+config['inj_nopb']
    inj_pb= sim_catalogs_path+config['inj_pb']
    rec_nopb= sim_catalogs_path+config['rec_nopb']
    rec_pb= sim_catalogs_path+config['rec_pb']
    _,header,w, BMAJ,_,_=correct_axes(residual_image)
    min_size=BMAJ*3600 #arcsec
    max_size=3*BMAJ*3600 
    survey_area=area**2*np.pi
    Svals = np.linspace(min_flux, max_flux, 10000)
    number_sources= tot_num_sources(Svals, survey_area)
    output_cat=path+catalogs
    #catalog_generation(pb_corr_image, name)
    #simulation(residual_image,pb_residual_image, min_size,max_size,sim_images_path,no_of_simulations, number_sources, Svals) # Perform simulation and injected catalogs
    #recovered_catalogs(sim_images_path,sim_catalogs_path,no_of_simulations) # create recovered catalogs
    #merg_catalogs(sim_catalogs_path)
    plots(merged_rec_nopb,merged_rec_pb,inj_nopb,inj_pb,rec_nopb,rec_pb,sim_plots_path)






