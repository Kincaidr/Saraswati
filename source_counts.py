import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def scale_flux(flux):
 counts_freq=1400
 Spectral_Index=-0.7
 data_freq=325
 flux=flux*(counts_freq/data_freq)**Spectral_Index
 return(flux)

def SCs_CLASS(S,  a0=3.192, a1=-0.223, a2=-0.846, a3=-0.261, a4=-0.024):
    a=np.array([a0, a1, a2, a3, a4])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
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

def norm(source_num,difflr,survey_area,centres):
    source_tot1=source_num/ (difflr)
    source_tot2=source_tot1/ ((survey_area)*(np.pi/180)**2)
    source_norm=source_tot2*(centres)**(2.5)
    return(centres, difflr, source_num, source_norm)

def resolution_bias_correction(resolution_bias_filename,centres):
    data = np.loadtxt(resolution_bias_filename)
    ratio=data[:,0]
    flux_extrap=data[:,1]
    f=interp1d(flux_extrap, ratio, bounds_error=False, fill_value="extrapolate")
    correction_factor=f(centres*1e3)
    return correction_factor    

def completeness_correction(output_filename, centres):
    data = np.loadtxt(output_filename)
    ratio=data[:,2]
    flux_extrap=data[:,0]*1e-3
    f=interp1d(flux_extrap, ratio, bounds_error=False, fill_value="extrapolate")
    correction_factor=f(centres)
    #correction_factor=np.append(correction_factor,np.ones(len(centres) - len(correction_factor)))
    corr=1/correction_factor
    return(corr)    

def false_detection_correction(output_filename):
    data = np.loadtxt(output_filename)
    corr=data[:,1]
    return(corr)    

def get_counts(real_cat, nbins):
    rec_cat = Table.read(real_cat)
    #mask=(rec_cat['S_Code'] =='S') #| (rec_cat['S_Code'] =='C') #| (rec_cat['S_Code'] =='M')
    flux=rec_cat['Total_flux']
    flux=flux*(counts_freq/data_freq)**Spectral_Index
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(0.03), num=nbins+1)
    #Range_x = np.percentile(flux, np.linspace(0, 80, nbins + 1))
    centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
    difflr=np.diff(Range_x)
    hist, x = np.histogram(flux, bins=Range_x)
    centres, deltaS, num, counts_norm=norm(hist,difflr,survey_area,centres)
    _,_,_,counts_err=norm(np.sqrt(hist),difflr,survey_area,centres)  
    return(centres, counts_norm,counts_err,deltaS, num)
    

def corrections(numbers,deltaS,centres,counts_norm,counts_err,source_count_literature,source_count_literature_DEEP,source_count_literature_SuperCLASS,completeness_corr,false_detect_corr, resolution_corr):
    # x_LoTSS=LoTSS_lit[:,0]
    # y_LoTSS=LoTSS_lit[:,1]
    # y_LoTSS_err_up=LoTSS_lit[:,2]
    # y_LoTSS_err_down=LoTSS_lit[:,3]
    SuperCLASS_lit=np.loadtxt(source_count_literature_SuperCLASS)
    x1_superCLASS=SuperCLASS_lit[:,0]
    x2_superCLASS=SuperCLASS_lit[:,1]
    xc_superCLASS=SuperCLASS_lit[:,2]
    N_superCLASS=SuperCLASS_lit[:,3]
    flux1=scale_flux(x1_superCLASS)*1e-3
    flux2=scale_flux(x2_superCLASS)*1e-3
    flux_xc=scale_flux(xc_superCLASS)*1e-3
    flux_xc=(flux1+flux2)/2
    difflr_superCLASS=flux2-flux1
    survey_area_CLASS=6.5
    centres_CLASS, deltaS_CLASS, num, counts_norm_CLASS=norm(N_superCLASS,difflr_superCLASS,survey_area_CLASS,flux_xc)
    DEEP_lit=np.loadtxt(source_count_literature_DEEP)
    x_DEEP=(10**DEEP_lit[:,0])
    y_DEEP=(10**DEEP_lit[:,1])*np.sqrt(x_DEEP)
    x_err_DEEP=DEEP_lit[:,2]
    y_err_DEEP=DEEP_lit[:,2]

    lit=np.loadtxt(source_count_literature)
    x_lit=(10**lit[:,0])*1e3
    y_lit=lit[:,1]
    x_err=lit[:,2]
    y_err=lit[:,3]
    var=lit[:,4]
    S=np.logspace(-2,3,1000) 
    plt.figure(figsize=(10, 7)) 
    plt.plot(S,10**SCs_Bondi(S),label='COSMOS 1.4GHz Bondi fitted (Bondi +2008)')
    #plt.plot(S,10**SCs_CLASS(S*1e-3),label='Super-CLASS Risely fitted (Risely +2016)')
    #plt.scatter(centres_CLASS*1e3,counts_norm_CLASS,label='SuperCLASS GMRT 325MHz (Risely +2016)',marker='o',color='red')
    plt.errorbar(centres[0]*1e3, counts_norm[0],yerr=counts_err[0],color='#ADD8E6' ,fmt='o',markersize=4,alpha=1,label=f'uncorrected A2631 MeerKAT 1.28GHz (This paper)')
    plt.errorbar(centres[1]*1e3, counts_norm[1],yerr=counts_err[1],color='#F08080' ,fmt='o',markersize=4,alpha=1,label=f'uncorrected Zwcl2341 MeerKAT 1.28GHz (This paper)')
    plt.errorbar(centres[0]*1e3, counts_norm[0]*completeness_corr[0],yerr=counts_err[0],xerr=deltaS[0],color='blue',label=f'A2631 MeerKAT 1.28GHz (This paper)' ,fmt='+',markersize=8,alpha=1)   
    plt.errorbar(centres[1]*1e3, counts_norm[1]*completeness_corr[1],yerr=counts_err[1],xerr=deltaS[1],color='red',label=f'Zwcl2341 MeerKAT 1.28GHz (This paper)' ,fmt='+',markersize=8,alpha=1) 
    #plt.errorbar(centres*1e3, counts_norm*completeness_corr,yerr=counts_err,color='orange',label=f'This paper corrected source counts' ,fmt='+')   
    plt.errorbar(x_lit,y_lit,yerr=[y_err+var,x_err-var], fmt='*',color='purple',label='de Zotti 1.4GHz compilation (de Zotti +2010)',alpha=0.3)
    plt.errorbar(x_DEEP*1e3,y_DEEP,yerr=[y_err_DEEP,x_err_DEEP], fmt='*',color='orange',label='DEEP2 MeerKAT 1.28GHz  (Mauch +2019)',alpha=0.7)
    #plt.errorbar(x_LoTSS,y_LoTSS,yerr=[y_LoTSS_err_up, y_LoTSS_err_down], fmt='*',color='orange',label='LoTSS 150MHz  (Mandal +2019)',alpha=0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(10**-3,10**4)
    plt.xlim(10**-2,10**2.5)
    plt.xlabel(r'$S_{1.4 GHz}$ [mJy ]',fontsize=18)
    plt.ylabel(r'$S^{2.5}dN/dS$ [sr$^{-1}$ Jy$^{1.5}$]',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=15, length=3, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=15, length=3, width=1)
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/Source_Counts_cone.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
    plt.close()

    with open("A2631_source_counts.txt", "w") as file:
        for c, d, n, counts, err, x, y, z in zip(centres[0],deltaS[0], numbers[0],counts_norm[0], counts_err[0], completeness_corr[0],resolution_corr[0], 1-false_detect_corr[0] ):
            file.write(f"{c*1e3} {d*1e3} {n} {counts} {err} {x} {y} {z}  \n")
    with open("Zwcl2341_source_counts.txt", "w") as file:
        for c, d, n, counts, err, x, y,z in zip(centres[1],deltaS[1], numbers[1],counts_norm[1], counts_err[1], completeness_corr[1],resolution_corr[1],1-false_detect_corr[1] ):
            file.write(f"{c*1e3} {d*1e3} {n} {counts} {err} {x} {y} {z}  \n")

if __name__ == "__main__":
    outname='/home/kincaid/Desktop/Saraswati_codes/plots/'
    source_count_literature='de_zotti.txt'
    source_count_literature_DEEP='MeerKAT_deep.txt'
    source_count_literature_SuperCLASS='SuperCLASS_source_counts.txt'
    survey_area=(45/60)*(45/60)
    nbins=20
    counts_freq=1.4
    data_freq=1.28  
    Spectral_Index=-0.7
    names=['A2631','Zwcl2341']
    centres_array=[]
    counts_norm_array=[]
    counts_err_array=[]
    completeness_corr_array=[]
    resolution_corr_array=[]
    false_detect_corr_array=[]
    deltaS_array=[]
    numbers_array=[]
    for name in names:
        resolution_bias_file=name+'_resolution_bias_correction.txt'
        completeness_file = name+"_visib_correction.txt"
        false_detection_file='false_detection_correction.txt'
        real_catalog =  '/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_srl.fits'   
        print( 'Real catalog is',real_catalog )
        centres,counts_norm,counts_err,deltaS, num=get_counts(real_catalog,nbins)
        numbers_array.append(num)
        deltaS_array.append(deltaS)
        centres_array.append(centres)
        counts_norm_array.append(counts_norm)
        counts_err_array.append(counts_err)
        completeness_corr= completeness_correction(completeness_file,centres)
        completeness_corr_array.append(completeness_corr)
        resolution_corr=resolution_bias_correction(resolution_bias_file, centres)
        resolution_corr_array.append(resolution_corr)
        false_detect_corr=false_detection_correction(false_detection_file)
        false_detect_corr_array.append(false_detect_corr)
    corrections(numbers_array,deltaS_array,centres_array,counts_norm_array,counts_err_array,source_count_literature, source_count_literature_DEEP,source_count_literature_SuperCLASS,completeness_corr_array, false_detect_corr_array, resolution_corr_array )  #calcaulte the ratio of recovered/injected in each bin
    
    
