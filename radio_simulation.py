#!/usr/bin/env python
# robert.kincaid@epfl.ch

#Create radio catalog from simulated radio sources with known gaussian noise
import os
import numpy as np
from astropy.table import Table
from astropy.io import fits
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
import subprocess

def false_Detections(catalog):
    cat= Table.read(catalog)
    num=len(cat)
    return(num)

def delete_simulation_contents(cluster_name):
    bash_script_path = "./delete_contents.sh"
    result = subprocess.run(
        [bash_script_path, cluster_name],
        text=True,          # Captures output as string
        capture_output=True  # Captures stdout and stderr
    )
    print(result.stdout)

def catalog_generation(fits_image, output_path, outname, res_image=True):
    outfile=output_path+outname+'_srl.fits'
    img = bdsf.process_image(fits_image, rms_box=(40,40),rms_box_bright=(20,20),adaptive_thresh=150,thresh_isl=4.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = False)
    img.write_catalog(outfile=outfile,format='fits', catalog_type='srl',clobber=True)
    print("Real catalog written")

    if res_image:
        img.export_image(outfile=output_path+outname+"_res_map.fits",clobber=True,img_type='gaus_resid')
        img.export_image(outfile=output_path+outname+"_rms_map.fits",clobber=True,img_type='rms')
        print("residual image written")
    else:
        print("residual image not written")

    return(outfile)


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


def place_sources_circ(num_sources, radius, min_distance, pix_size, centre_x, centre_y):
    sources = []
    attempts = 0
    max_attempts = num_sources * 10  # Limit attempts to avoid infinite loop
    
    while len(sources) < num_sources and attempts < max_attempts:
        # Generate random angle and distance for a point within the circle
        angle = np.random.uniform(0, 2 * np.pi)
        distance = np.sqrt(np.random.uniform(0, 1)) * radius
        xc = int(centre_x + distance * np.cos(angle))
        yc = int(centre_y + distance * np.sin(angle))
        new_source = (xc, yc)
        if not is_too_close(new_source, sources, min_distance, pix_size):
            sources.append(new_source)
        attempts += 1
    if len(sources) < num_sources:
        print(f"Warning: Only placed {len(sources)} sources out of {num_sources} after {max_attempts} attempts.")
    return sources


def place_sources_box(num_sources, size, min_distance, pix_size, centre_x, centre_y):
    sources = []
    attempts = 0
    max_attempts = num_sources * 10  # Limit attempts to avoid infinite loop
    size=size*60
    width= int(size / pix_size)  # Convert width from arcsec to pixels
    height= int(size / pix_size)  # Convert height from arcsec to pixels

    half_width = width // 2
    half_height = height // 2

    while len(sources) < num_sources and attempts < max_attempts:
        x = int(np.random.uniform(centre_x - half_width, centre_x + half_width))
        y = int(np.random.uniform(centre_y - half_height, centre_y + half_height))
        new_source = (x, y) 
        if not is_too_close(new_source, sources, min_distance, pix_size):
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

def cross_match(cat1,cat2,seperation,RA_1,DEC_1,RA_2,DEC_2):
        DEC_inj=cat1[DEC_1]
        DEC_rec=cat2[DEC_2]
        RA_inj=cat1[RA_1]
        RA_rec=cat2[RA_2]
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

def image_properties(fits_image):
    data_hdu = fits.open(fits_image)[0]
    data_data = data_hdu.data
    data_header = data_hdu.header
    BMAJ=data_header['BMAJ']
    BMIN=data_header['BMIN']
    BPA=data_header['BPA']
    centre_x=data_header['CRVAL1']
    centre_y=data_header['CRVAL2']
    new_data = data_data.squeeze()
    new_data=new_data.T
    w = WCS(fits_image, relax=True)
    return(new_data,data_header,w, BMAJ,BMIN,BPA, centre_x, centre_y)

def convolve_gaussian( BMAJ, pix_size):
    fwhm_x=np.random.choice(np.random.uniform(BMAJ, 2*BMAJ,1))
    fwhm_y=np.random.choice(np.random.uniform(0.5,1,1))*fwhm_x
    sigma_x = fwhm_x / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    sigma_y = fwhm_y / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    return(fwhm_x, fwhm_y, sigma_x, sigma_y)

def convolve_gaussian_unresolved( BMAJ, pix_size):
    fwhm_x=BMAJ
    fwhm_y=fwhm_x
    sigma_x = fwhm_x / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    sigma_y = fwhm_y / (2.0 * np.sqrt(2.0 * np.log(2.0))*pix_size)
    return(fwhm_x, fwhm_y, sigma_x, sigma_y)

def gaussianv1(A, sigma_x,sigma_y, theta):
    sigma = 10
    size1 = int(np.round((sigma_x))*sigma)
    size2 = int(np.round((sigma_y))*sigma)
    x = np.arange(0, size1, dtype=np.float64)[:, None]
    y = np.arange(0, size2, dtype=np.float64)[None, :]
    gaussian_model = Gaussian2D(amplitude=A, x_mean=size1/2, y_mean=size2/2, x_stddev=sigma_x, y_stddev=sigma_y, theta=theta)
    gaussian_source = gaussian_model(x, y)
    return gaussian_source

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def distance(x1,x2,y1,y2):
    dist=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dist

def get_source_fluxes(min_flux,max_flux, numpoints):
    Svals = np.linspace(min_flux, max_flux, 10000)
    dNdSS25 = 10 ** (SCs_Bondi(Svals))
    flux_prob = dNdSS25 * (Svals/1000)**(-2.5)
    flux_prob_norm=flux_prob/np.sum(flux_prob)
    CDF=np.cumsum(flux_prob_norm)
    inverse_CDF = interpolate.interp1d(CDF, Svals, bounds_error=False, fill_value=(Svals[0], Svals[-1]))
    x=np.random.uniform(0,1,numpoints)
    flux_samples=inverse_CDF(x)
    return flux_samples


def source_properties_resolved(xc, yc, flux, BMAJ,BMIN, BPA, min_size, max_size, pix_size, w, padding):
        fwhm_x, fwhm_y, sigma_x, sigma_y=convolve_gaussian( BMAJ, pix_size)
        area_source= fwhm_x* fwhm_y
        area_beam=BMAJ*BMIN
        A=flux*(area_beam/area_source)    
        BPA=np.random.randint(0,180)
        gaussian=gaussianv1(A,sigma_x, sigma_y, BPA)
        x1=int(xc-(gaussian.shape[0]/2))
        x2=int(xc+(gaussian.shape[0]/2))
        y1=int(yc-(gaussian.shape[1]/2))
        y2=int(yc+(gaussian.shape[1]/2))
        xc, yc=xc-padding/2, yc-padding/2   
        aux_deg = w.pixel_to_world(yc, xc, 1,StokesCoord(1))[0]
        xc_deg, yc_deg = aux_deg.ra.value, aux_deg.dec.value  
        return(x1,x2,y1,y2,xc_deg, yc_deg, gaussian,fwhm_x, fwhm_y,A, flux)


def source_properties_unresolved(xc, yc, flux, BMAJ,BMIN, BPA, pix_size, w, padding):
        BMAJ=BMIN
        fwhm_x, fwhm_y, sigma_x, sigma_y=convolve_gaussian_unresolved( BMAJ, pix_size)
        A=flux
        gaussian=gaussianv1(A,sigma_x, sigma_y, BPA)
        x1=int(xc-(gaussian.shape[0]/2))
        x2=int(xc+(gaussian.shape[0]/2))
        y1=int(yc-(gaussian.shape[1]/2))
        y2=int(yc+(gaussian.shape[1]/2))
        xc, yc=xc-padding/2, yc-padding/2   
        aux_deg = w.pixel_to_world(yc, xc, 1,StokesCoord(1))[0]
        xc_deg, yc_deg = aux_deg.ra.value, aux_deg.dec.value  
        return(x1,x2,y1,y2,xc_deg, yc_deg, gaussian,fwhm_x, fwhm_y, A, flux)

def injected_catalogs(min_size,max_size,fluxes,sim_images_path, number_sources, source_placement, data, output) :
    real_noise=data
    Data_real=np.zeros(padded_image.shape)
    RA=np.zeros(number_sources)
    DEC=np.zeros(number_sources)
    RA_pix=np.zeros(number_sources)
    DEC_pix=np.zeros(number_sources)
    Peak=np.zeros(number_sources)
    Maj_conv=np.zeros(number_sources)
    Min_conv=np.zeros(number_sources)
    Int_size=np.zeros(number_sources)
    Flux_total=np.zeros(number_sources)
    Dist=np.zeros(number_sources)
    number_sources_resolved=int(0.4*number_sources)

    for i in range(number_sources_resolved): 
        xc = source_placement[i][0]
        yc = source_placement[i][1]                                                                                      
        x1,x2,y1,y2,xc_deg, yc_deg, gaussian,maj_conv, min_conv, A, Total_flux=source_properties_resolved(xc, yc, fluxes[i]*10**-3, BMAJ,BMIN,BPA, min_size, max_size, pix_size,w, padding)
        Data_real[x1:x2,y1:y2] +=gaussian
        RA[i]=xc_deg-360
        DEC[i]=yc_deg
        RA_pix[i]=xc
        DEC_pix[i]=yc
        Peak[i]=A
        Dist[i]=distance(RA[i],RA_centre,DEC[i],DEC_centre)
        Maj_conv[i]=maj_conv
        Min_conv[i]=min_conv
        Flux_total[i]=Total_flux
        print('Source ' + str(i) + ' written')

    for i in range(number_sources_resolved, number_sources): 
        xc = source_placement[i][0]
        yc = source_placement[i][1]                                                                                      
        x1,x2,y1,y2,xc_deg, yc_deg, gaussian,maj_conv, min_conv, A, Total_flux=source_properties_unresolved(xc, yc, fluxes[i]*10**-3, BMAJ,BMIN,BPA, pix_size,w, padding)
        Data_real[x1:x2,y1:y2] +=gaussian
        RA[i]=xc_deg-360
        DEC[i]=yc_deg
        RA_pix[i]=xc
        DEC_pix[i]=yc
        Peak[i]=A
        Dist[i]=distance(RA[i],RA_centre,DEC[i],DEC_centre)
        Maj_conv[i]=maj_conv
        Min_conv[i]=min_conv
        Flux_total[i]=Total_flux
        print('Source ' + str(i) + ' written')

    final_image=Data_real[paddiv2:-paddiv2,paddiv2:-paddiv2]
    sim_image_test=sim_images_path+'simulated_image_test_'+output+'.fits'
    sim_image_real=sim_images_path+'simulated_image_real_'+output+'.fits'
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
    col7 = fits.Column(name='Int_size', format='D',array=Int_size)
    col8 = fits.Column(name='Maj_conv', format='D',array=Maj_conv)
    col9 = fits.Column(name='Min_conv', format='D',array=Min_conv)
    col10 = fits.Column(name='Dist', format='D',array=Dist)
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10])
    inj_cat_real=sim_catalogs_path+'injected_cat_'+ output+'.fits'
    hdu.writeto(inj_cat_real,overwrite=True)
    print(inj_cat_real+ ' ' +'written')
    return inj_cat_real, sim_image_real

def recovered_catalogs(sim_image, output):
    output_name=output
    output_path=sim_catalogs_path
    outfile=catalog_generation(sim_image, output_path, output_name, res_image=False)
    print(f'Real catalog written: {str(output_name)}')
    return outfile
 

def update_table_v2(filename, count, centre, sim_no, nbins):
    try:
        with open(filename, "r") as file:
            lines = file.readlines()
            bin_centres = lines[0].strip().split("\t")  # First row: bin centres
            counts_existing = [line.strip().split("\t") for line in lines[1:]]  # Other rows

    except FileNotFoundError:
        bin_centres = []
        counts_existing = []

    if f"{centre:.5f}" not in bin_centres:
        bin_centres.append(f"{centre:.5f}")

    column_index = bin_centres.index(f"{centre:.5f}")
    while len(counts_existing) <= sim_no:
        counts_existing.append(["0"] * nbins)  # Add enough rows for sim_no

    counts_existing[sim_no][column_index] = str(count)   
    with open(filename, "w") as file:
        file.write("\t".join(bin_centres) + "\n")  # First row: bin centres
        for row in counts_existing:
            file.write("\t".join(row) + "\n")  # Other rows for simulations
    print(f"Updated {filename} with Bin {centre:.5f}, Count {count}, Simulation {sim_no}.")

def ratio(rec_cat, centre, sim_no, nbins, false_detections, rec_counts_file):
    recovered_cat = Table.read(rec_cat)
    rec_count = len(recovered_cat['Total_flux'])-false_detections
    update_table_v2( rec_counts_file, rec_count, centre, sim_no, nbins)


def simulation(min_size,max_size,sim_images_path, number_sources, Range_x,data, rec_counts_file):
    if os.path.exists(rec_counts_file):
        os.remove(rec_counts_file)
        print(f"Removed {rec_counts_file}")
    else:
        print(f"{rec_counts_file} does not exist.")
    for sim_no in range(no_of_simulations):
        source_placement=place_sources(number_sources, xx, yy, min_distance=20, pix_size=pix_size)
        print("Simulation_number is", sim_no)
        for bin_no in range(len(Range_x)-1):
            left_bin=Range_x[bin_no]
            right_bin=Range_x[bin_no+1]
            centre=(right_bin+left_bin)/2
            flux_samples=get_source_fluxes(left_bin,right_bin,number_sources)
            output='sim_'+str(sim_no)+'_bin_'+str(bin_no)
            print("run for "+ output)
            injected_catalog, simulated_image=injected_catalogs(min_size,max_size,flux_samples,sim_images_path, number_sources,source_placement, data, output) # Perform simulation and injected catalogs
            recovered_catalog=recovered_catalogs(simulated_image, output) 
            ratio(recovered_catalog, centre, sim_no, nbins, false_detections, rec_counts_file)
##
if __name__ == "__main__":
    name='Zwcl2341'
    min_size = 1  # Minimum source size in simulation (arcsec)
    max_size = 30 # Maximum source size in simulation (arcsec)
    number_sources=500 # total number of sources simulated on each image
    no_of_simulations=1 # total number of simulations
    flux_min=0.03#minimum flux of sources in simualtion (mJy)
    flux_max=10 #maximum flux of sources in simualtion (mJy)
    nbins=30 #Number of bins you want to find the detected fraction
    residual_image = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'+name+'_cut_res_map.fits'
    residual_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_cut_srl.fits' #residual image for which we place sources on
    sim_images_path = '/media/kincaid/LaCie/images/' #all simulated images go here
    sim_catalogs_path = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/simulation/sim_catalogs/' #all simulated catalogs go here
    Range_x = 10**np.linspace(start=np.log10(flux_min), stop=np.log10(flux_max), num=nbins+1)
    padding=1000
    data,header,w, BMAJ,BMIN,BPA, RA_centre, DEC_centre=image_properties(residual_image) #extract fits data, header 
    original_image = np.zeros((header['NAXIS1'], header['NAXIS2']))
    padded_image= np.zeros((original_image.shape[0]+padding, original_image.shape[1]+padding))
    noise=1e-9
    #radius=1680
    mean,stddev=0,noise
    uniform_noise = np.random.normal(mean, stddev, original_image.shape)
    pix_size=header['CDELT2']*3600
    BMAJ=BMAJ*3600
    BMIN=BMIN*3600
    rec_counts_file = name+"_recovered_counts_table_cut.txt"
    #This section is for circular cutout image, need to uncomment sources_placement_circle in simualtion function to use
    paddiv2=int(padding/2)
    xx=(paddiv2,paddiv2+original_image.shape[0])
    yy=(paddiv2,paddiv2+original_image.shape[1])
    centre_x = (paddiv2 + original_image.shape[0] / 2) 
    centre_y = (paddiv2 + original_image.shape[1] / 2) 
    #radius = 1650 # circular cutout size radius in pixels
    ##########################################################

    #res_cat=catalog_generation(residual_image, sim_catalogs_path, outname=name+'_res', res_image=False) #Run source finder on residual image to find false detections
    false_detections=0#106 205q#false_Detections(res_cat) #the number of false detections, can also specify manually
    #print("False detections are", false_detections) #   specifiy correct false detections here before simulation
    simulation(min_size,max_size,sim_images_path, number_sources, Range_x, data, rec_counts_file)


  



    


