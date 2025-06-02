import numpy as np
from astropy.table import Table, vstack
from astropy.io import fits
import matplotlib.pyplot as plt
import json

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals


def norm(source_num,difflr,survey_area,centres):
    source_tot1=source_num/ (difflr)
    source_tot2=source_tot1/ ((survey_area)*(np.pi/180)**2)
    source_norm=source_tot2*centres**(2.5)
    return(source_norm)

def source_counts_v2(flux,survey_area,counts_freq,data_freq, Spectral_Index,nbins,corr=None,Range_x=None):
        
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


def ratio_source_count(flux, nbins, equal_source_size=False, equal_bin_size=False, Range_x= None):

        if  equal_source_size:
                pers = np.linspace(0, 100, nbins+1)
                Range_x = np.percentile(flux, pers)
        elif equal_bin_size:
                Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(flux.max())+1e-8, num=nbins+1)
        else:
              Range_x=Range_x
        centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
        hist, _ = np.histogram(flux, bins=Range_x)
        source_tot=hist
        return(Range_x, source_tot, centres)

def corrections(merged_inj,merged_rec,sigma,output_cat ,output_img):
   
    rec_cat=merged_rec
    inj_cat_real = Table.read(merged_inj)
    rec_cat_real = Table.read(rec_cat)
    flux_inj=inj_cat_real['Total_flux_inj']
    flux_rec=rec_cat_real['Total_flux'][rec_cat_real['Total_flux'] < 0.01]
    nbins=50
    
    #Range_x=np.array([2e-05,5e-05,6e-05,6.5e-05,7e-05,7.5e-05,8e-05,8.5e-05,9e-05,9.5e-05,1e-04,1.1e-04,1.2e-04,1.3e-04,1.4e-04,1.5e-04,1.6e-04,1.8e-04,1.9e-04,2e-04,2.2e-04,2.4e-04,2.6e-04,3e-04,4e-04,4.5e-04,5e-04,8e-04,1e-03,2e-03,1e-02])    
    Range_x,counts_rec,_=ratio_source_count(flux_rec, nbins=nbins, equal_source_size=True, equal_bin_size=False, Range_x= None)
    _,counts_inj, bin_centres=ratio_source_count(flux_inj, nbins=nbins, equal_source_size=False, equal_bin_size=False, Range_x= Range_x)

    breakpoint()
    ratio=counts_rec/counts_inj
    error=ratio*np.sqrt((np.sqrt(counts_inj)/counts_inj)**2+(np.sqrt(counts_rec)/counts_rec)**2)
    print(ratio)
    plt.errorbar(bin_centres*10**3,ratio,yerr=error)
    plt.xlabel('Log Flux S [mJy]',size=12)
    plt.ylabel('Recovered/Injected',size=12)
    plt.axvline(x=5*sigma, color='red', linestyle='--', label=r'detection threshold $5 \sigma$')
    plt.axhline(y=1,color='green')
    plt.legend()
    plt.xscale('log')
    plt.savefig(output_img+'ratio_plot.png')
    plt.show()
    plt.close()
    col1 = fits.Column(name='ratio', format='D',array=ratio)
    col2 = fits.Column(name='error', format='D',array=error)
    hdu = fits.BinTableHDU.from_columns([col1, col2])
    correction_cat=output_cat+'incompleteness_catalog.fits'
    hdu.writeto(correction_cat,overwrite=True)
    print(correction_cat+ ' written')

def source_count_corrections(incompleteness_cat,flux_corr_catalog,survey_area,output_cat,output_img ):

    corrections_cat=incompleteness_cat
    flux=Table.read(flux_corr_catalog)['Total_flux']  
    corrections_cat = Table.read(corrections_cat)
    nbins=len(corrections_cat)
    ratio=np.array(corrections_cat['ratio'])
    #ratio[0]=ratio[0]*10**2
    correction=[1/x if x != 0 else float('inf') for x in ratio]
    Range_x=np.array([2e-05,5e-05,6e-05,6.5e-05,7e-05,7.5e-05,8e-05,8.5e-05,9e-05,9.5e-05,1e-04,1.1e-04,1.2e-04,1.3e-04,1.4e-04,1.5e-04,1.6e-04,1.8e-04,1.9e-04,2e-04,2.2e-04,2.4e-04,2.6e-04,3e-04,4e-04,4.5e-04,5e-04,8e-04,1e-03,2e-03,1e-02])
    _,_,counts,counts_err,bin_centre=source_counts_v2(flux=flux,survey_area=survey_area,data_freq=1.28,counts_freq=1.4, Spectral_Index=-0.7,nbins=nbins, Range_x=None)
    _,_,counts_corr,counts_corr_err,_=source_counts_v2(flux=flux,survey_area=survey_area,data_freq=1.28,counts_freq=1.4, Spectral_Index=-0.7,nbins=nbins,corr=correction,Range_x=None)
    plt.errorbar(bin_centre*10**3, counts,yerr=counts_err,color='red',label=f'Uncorrected source counts' ,fmt='o')
    plt.errorbar(bin_centre*10**3, counts_corr,yerr=counts_corr_err,color='orange',label=f'Corrected source counts' ,fmt='+')   
    data_to_save = np.column_stack((bin_centre, counts, counts_err, counts_corr, counts_corr_err))
    np.savetxt(output_cat+'source_counts.txt', data_to_save, header='bin_centre counts counts_err counts_corr counts_corr_err', fmt='%f')

    S=np.logspace(-1.5,3,10001)
    plt.plot(S,10 ** (SCs_Bondi(S)),label='Bondi 2008 1.4 GHz')
    
    #plt.plot(S,10 ** (SCs_Risely(S)),label='Risely 2016 1.4 GHz')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(10**-3,10**3)
    plt.xlabel('Flux S [mJy]',size=12)
    plt.ylabel('S^(2.5) * dN/dS',size=12)
    plt.title('Source Counts')
    plt.legend()
    plt.savefig(output_img+'Source_Counts.png')
    plt.show()
    plt.close()

def read_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)
    
if __name__ == "__main__":

    config_path = 'A2631.json'
    config = read_config(config_path)
    name=config['name']
    path=config['path']
    output_catalogs=config['catalogs']
    output_plots=config['plots']
    merged_inj= output_catalogs+config['merged_inj']
    merged_rec= output_catalogs+config['merged_rec']
    merged_rec_corr= config['merged_rec_corr']
    incompleteness_cat= config["incompleteness_cat"]
    nbins = config['nbins']
    sigma = config['sigma']
    area=config["area"]
    survey_area=area**2*np.pi
    real_catalog =  config['real_catalog']
    corrections(merged_inj,merged_rec,sigma, output_catalogs, output_plots)  #calcaulte the ratio of recovered/injected in each bin
    source_count_corrections(incompleteness_cat,real_catalog,survey_area,output_catalogs, output_plots ) # Apply correction to the real source counts
    