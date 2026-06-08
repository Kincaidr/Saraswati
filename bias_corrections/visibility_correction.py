from astropy.io import fits
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def visib_plot(real_cat, fits_image,name):
    cat = Table.read(real_cat)
    flux=cat['Total_flux']*1e3
    hdu = fits.open(fits_image)[0]
    image_data = hdu.data*1e3

    if len(image_data.shape) > 2:
        image_data = image_data[0, 0,:, :] 
   
    mask=np.isnan(image_data)
    image_data=image_data[~mask]
    header = fits.getheader(fits_image)
    flux_min=0.01
    flux_max=10
    irange = 10**np.linspace(start=np.log10(flux_min), stop=np.log10(flux_max), num=nbins+1)
    centres = (irange[0:-1] + irange[1:]) / 2.0
    area_i = np.ones(len(irange))
    #pixels_tot = np.sum(image_data < irange[-1]/detect_thresh)
    pixels_tot = np.sum(image_data < irange[-1])
    Area_tot =  pixels_tot
    rms_lim_array=[]

    for i in range(len(irange)-1):
        rmslim = np.sqrt(irange[i]*irange[i+1])/(5)
        #rmslim = np.sqrt(irange[i]*irange[i+1])
        pixels = np.sum(image_data < rmslim)
        rms_lim_array.append(rmslim)
    
        if pixels == 0:
            area_i[i] = np.nan
        else:
            Area = pixels
            area_i[i] = Area/Area_tot
    
    correcion_file=name+"_corrected_area.txt"
    with open(correcion_file, "w") as file:
        for value, flux in zip(area_i,  centres):
            file.write(f"{1/value} {flux}  \n")
    return(correcion_file)

    
def completeness_file(output_filename,inj_sources, no_simulations ):
    data = np.loadtxt(output_filename )
    bin_centres = data[0, :]
    old_counts = data[1:, :] 
    new_counts = np.zeros(old_counts.shape)
    for i in range(len(old_counts[:,0])):
        for j in range(len(old_counts[0,:])):
            new_counts[i][j] = old_counts[i][j]-(false_sources)

    rec_counts = np.sum(new_counts[0:, :], axis=0)
    corr =(np.max(rec_counts) -(inj_sources*no_simulations)) # Ensure no negative counts
    rec_counts = rec_counts - corr
    print("Bin Centres:", bin_centres)
    print("Summed Recovered Counts:", rec_counts)
    inj_sources_tot=np.full(nbins,inj_sources*no_simulations)
    ratio = rec_counts /inj_sources_tot
    error=np.sqrt( (np.sqrt(inj_sources)/inj_sources)**2+ (np.sqrt(rec_counts)/rec_counts)**2)

    output_file = name+"_output_table.txt"
    with open(output_file, "w") as file:
        for r, rat, err in zip(bin_centres,ratio, error ):
            file.write(f"{r} {rat} {err}\n")
    print(f"Ratio array has been written to {output_file}")
    return(output_file)


# def get_data(visibility_correction_file,  completeness):    
#     comp_file=np.loadtxt(completeness)
#     visib_file=np.loadtxt(visibility_correction_file)
#     visib_scale=visib_file[:,0]
#     visib_bin=visib_file[:,1]
#     comp=comp_file[:,1]
#     comp_bin=comp_file[:,0]
#     comp_err=comp_file[:,2]
#     f=interp1d(visib_bin, visib_scale, bounds_error=False, fill_value="extrapolate")
#     visib_new=f(comp_bin)
#     comp_err=comp_err
#     comp_corr=comp#*visib_new
#     f=interp1d(comp_bin, comp_corr, bounds_error=False, fill_value="extrapolate")
#     comp_corr_new=f(comp_bin)
    
#     with open(name+'/'+'visib_correction_cut.txt', "w") as file:
#         for bin, co, co_corr, err in zip(comp_bin,comp,comp_corr_new, comp_err ):
#             file.write(f"{bin} {co} {co_corr} {err}\n")

def get_data( completeness):    
    comp_file=np.loadtxt(completeness)
    comp_bin=comp_file[:,0]
    comp=comp_file[:,1]
    comp_err=comp_file[:,2]
    f=interp1d(comp_bin, comp, bounds_error=False, fill_value="extrapolate")
    
    with open(name+'/'+'visib_correction_cut.txt', "w") as file:
        for bin, co, err in zip(comp_bin,comp, comp_err ):
            file.write(f"{bin} {co} {err}\n")

if "__name__":
    inj_sources=500
    no_simulations=20
    nbins=30
    name='Zwcl2341'
    false_sources=205 #106
    fits_image = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'+name+'_cut_rms_map.fits'   
    real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_srl_flux_corr.fits'
    completeness_filename = name+'/'+"recovered_counts_table_cut.txt"
    #visibility_correction_file=visib_plot(real_cat, fits_image, name)
    completeness=completeness_file(completeness_filename,inj_sources, no_simulations )
    get_data( completeness)




