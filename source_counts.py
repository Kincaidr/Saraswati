import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Massardi_2010_plot import Massardi_counts
from Mancuso_source_counts import Mancuso_counts
from SEMPER_source_counts import SEMPER_SFG_AGN_counts
from TRECS_source_counts import TRECS_counts
import pickle

def scale_flux(flux,data_freq,counts_freq=1400, Spectral_Index=-0.7):
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
    with open(resolution_bias_filename, 'rb') as f:
        func = pickle.load(f)
    corr=func(centres*1e3)
    breakpoint()
    return corr   

def completeness_correction(output_filename, centres):
    data = np.loadtxt(output_filename)
    flux_extrap=data[:,0]
    ratio=data[:,1]
    f=interp1d(flux_extrap, ratio, bounds_error=False, fill_value="extrapolate")
    correction_factor=f(centres*1e3)
    #correction_factor=np.append(correction_factor,np.ones(len(centres) - len(correction_factor)))
    corr=1/correction_factor
    corr=np.abs(corr)
    return(corr)    

def get_counts(real_cat, nbins):
    rec_cat = Table.read(real_cat)
    #mask=(rec_cat['S_Code'] =='S') #| (rec_cat['S_Code'] =='C') #| (rec_cat['S_Code'] =='M')
    flux=rec_cat['Total_flux']
    flux=flux*(counts_freq/data_freq)**Spectral_Index
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(0.2), num=nbins+1)
    #Range_x = np.percentile(flux, np.linspace(0, 80, nbins + 1))
    centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
    difflr=np.diff(Range_x)
    hist, x = np.histogram(flux, bins=Range_x)
    centres, deltaS, num, counts_norm=norm(hist,difflr,survey_area,centres)
    _,_,_,counts_err=norm(np.sqrt(hist),difflr,survey_area,centres)  
    return(centres, counts_norm,counts_err,deltaS, num)
    

def corrections(numbers,deltaS,centres,counts_norm,counts_err,source_count_literature,source_count_literature_DEEP,source_count_literature_SuperCLASS,completeness_corr, resolution_corr):
    # x_LoTSS=LoTSS_lit[:,0]
    # y_LoTSS=LoTSS_lit[:,1]
    # y_LoTSS_err_up=LoTSS_lit[:,2]
    # y_LoTSS_err_down=LoTSS_lit[:,3]
    SuperCLASS_lit=np.loadtxt(source_count_literature_SuperCLASS)
    MeerKAT_COSMOS=np.loadtxt(MeerKAT_COSMOS_lit)
    MeerKAT_XMM_LSS=np.loadtxt(MeerKAT_XMM_LSS_lit)
    x1_superCLASS=SuperCLASS_lit[:,0]
    x2_superCLASS=SuperCLASS_lit[:,1]
    xc_superCLASS=SuperCLASS_lit[:,2]
    N_superCLASS=SuperCLASS_lit[:,4]
    # flux1=scale_flux(x1_superCLASS)*1e-3
    # flux2=scale_flux(x2_superCLASS)*1e-3
    # flux_xc=scale_flux(xc_superCLASS)*1e-3
    # flux_xc=(flux1+flux2)/2
    # difflr_superCLASS=flux2-flux1
    # survey_area_CLASS=6.5
    #centres_CLASS, deltaS_CLASS, num, counts_norm_CLASS=norm(N_superCLASS,difflr_superCLASS,survey_area_CLASS,flux_xc)
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
    S_M,M_counts=Mancuso_counts()
    Semper_M,Semper_counts=  SEMPER_SFG_AGN_counts()
    TRECS_M,TRECS_count=  TRECS_counts()
    S=np.logspace(-2,3,1000) 
    fig, axs = plt.subplots(2, 1, figsize=(9, 12), sharex=True)

    # fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharey=True)

    # for i, name in enumerate(names):
    #     axs[i].set_xscale('log')
    #     axs[i].set_yscale('log')
    #     axs[i].set_ylim(1e-3, 1e4)
    #     axs[i].set_xlim(1e-2, 10**2.5)
    #     axs[i].set_xlabel(r'$S_{1.4 GHz}$ [mJy]', fontsize=16)
    #     axs[i].set_ylabel(r'$S^{2.5}dN/dS$ [sr$^{-1}$ Jy$^{1.5}$]', fontsize=16)
    #     axs[i].legend()
    #     axs[i].axvline(5*16e-3,linestyle='--',color='blue',alpha=1)
    #     axs[i].legend(loc='lower right')
    #     axs[i].tick_params(axis='both', which='major', labelsize=12)
    #     axs[i].tick_params(axis='both', which='minor', labelsize=12)
    #     axs[i].grid(True, which="both", ls="--", alpha=0.5)
    #     axs[i].grid(True, which="both", ls="--", alpha=0.5)
    #     axs[i].errorbar(centres[0]*1e3, counts_norm[0], yerr=counts_err[0], 
    #                 color='#ADD8E6', fmt='o', markersize=4, alpha=1,
    #                 label='MOSS2 A2631 uncorrected (This paper)')
    #     axs[i].errorbar(centres[0]*1e3, counts_norm[0]*completeness_corr[0], 
    #                 yerr=counts_err[0], xerr=deltaS[0],
    #                 color='blue', fmt='+', markersize=10, alpha=1,
    #                 label='MOSS2 A2631 corrected (This paper)')

    axs[0].grid(True, which="both", ls="--", alpha=0.5)
    axs[1].grid(True, which="both", ls="--", alpha=0.5)
    axs[0].errorbar(centres[0]*1e3, counts_norm[0], yerr=counts_err[0], 
                    color='#ADD8E6', fmt='o', markersize=4, alpha=1,
                    label='MOSS2 A2631 uncorrected (This paper)')

    axs[0].errorbar(centres[0]*1e3, counts_norm[0]*completeness_corr[0], 
                    yerr=counts_err[0], xerr=deltaS[0],
                    color='blue', fmt='+', markersize=16, alpha=1,
                    label='MOSS2 A2631 corrected (This paper)')
    axs[0].plot(S_M,M_counts,label='Mancuso 2017 model (Mancuso +2017)',color="#F700CD",alpha=0.7)
    axs[0].plot(Semper_M,Semper_counts,label='SEMPER 1.4GHz model (Giulietti +2025)',color='green',alpha=0.7)
    axs[0].plot(TRECS_M,TRECS_count,label='T-RECS model (Bonaldi +2023)',color='orange',alpha=0.7)
    axs[0].plot(S,10**SCs_Bondi(S),label='COSMOS 1.4GHz Bondi fitted (Bondi +2008)',alpha=0.7)
    axs[0].errorbar(x_lit,y_lit,yerr=[y_err+var,x_err-var], fmt='*',color='grey',label='de Zotti 1.4GHz compilation (de Zotti +2010)',alpha=0.2)
    axs[0].errorbar(x_DEEP*1e3,y_DEEP,yerr=[y_err_DEEP,x_err_DEEP], fmt='x',color='black',label='DEEP2 1.28GHz (Mauch +2019)',alpha=0.5,markersize=8)
    axs[0].errorbar(MeerKAT_XMM_LSS[:,0]*1e-3,MeerKAT_XMM_LSS[:,1],color="#BBBE16",yerr=[MeerKAT_XMM_LSS[:,2],MeerKAT_XMM_LSS[:,3]],markersize=6,alpha=0.4,label='MIGHTEE XMM LSS 1.4GHz (Hale +2022)' ,fmt='d')
    axs[0].errorbar(MeerKAT_COSMOS[:,0]*1e-3,MeerKAT_COSMOS[:,1],color='green',yerr=[MeerKAT_COSMOS[:,2],MeerKAT_COSMOS[:,3]],markersize=6,alpha=0.4,label='MIGHTEE COSMOS 1.4GHz (Hale +2022)' ,fmt='^')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylim(1e-3, 1e4)
    axs[0].set_xlim(1e-2, 10**2.5)
    axs[0].set_xlabel(r'$S_{1.4 GHz}$ [mJy]', fontsize=16)
    axs[0].set_ylabel(r'$S^{2.5}dN/dS$ [sr$^{-1}$ Jy$^{1.5}$]', fontsize=16)
    axs[0].legend()
    axs[0].axvline(5*16e-3,linestyle='--',color='blue',alpha=1)
    axs[0].legend(loc='lower right')
    axs[0].set_title('A2631', fontsize=18)
    axs[0].tick_params(axis='both', which='major', labelsize=12)
    axs[0].tick_params(axis='both', which='minor', labelsize=12)

    # ---- Plot 2: Zwcl2341 ----
    axs[1].errorbar(centres[1]*1e3, counts_norm[1], yerr=counts_err[1], 
                    color='#DA6969', fmt='o', markersize=4, alpha=1,
                    label='MOSS2 Zwcl2341 uncorrected (This paper)')
    axs[1].errorbar(centres[1]*1e3, counts_norm[1]*completeness_corr[1], 
                    yerr=counts_err[1]*resolution_corr[1], xerr=deltaS[1],
                    color='red', fmt='+', markersize=16, alpha=1,
                    label='MOSS2 Zwcl2341 corrected (This paper)')
    axs[1].plot(S_M,M_counts,label='Mancuso model (Mancuso +2017)',color="#F700CD",alpha=0.7)
    axs[1].plot(Semper_M,Semper_counts,label='SEMPER model (SEMPER +2025)',color='green',alpha=0.7)
    axs[1].plot(S,10**SCs_Bondi(S),label='COSMOS 1.4GHz (Bondi +2008)',alpha=0.7)
    axs[1].plot(TRECS_M,TRECS_count,label='T-RECS model (Bonaldi +2023)',color='orange')
    axs[1].errorbar(x_lit,y_lit,yerr=[y_err+var,x_err-var], fmt='*',color='grey',label='de Zotti 1.4GHz compilation (de Zotti +2010)',alpha=0.2)
    axs[1].errorbar(x_DEEP*1e3,y_DEEP,yerr=[y_err_DEEP,x_err_DEEP], fmt='x',color='black',label='DEEP2 1.28GHz (Mauch +2019)',alpha=0.5,markersize=8)
    axs[1].errorbar(MeerKAT_XMM_LSS[:,0]*1e-3,MeerKAT_XMM_LSS[:,1],color="#BBBE16",yerr=[MeerKAT_XMM_LSS[:,2],MeerKAT_XMM_LSS[:,3]],markersize=6,alpha=0.4,label='MIGHTEE XMM LSS 1.4GHz (Hale +2022)' ,fmt='d')
    axs[1].errorbar(MeerKAT_COSMOS[:,0]*1e-3,MeerKAT_COSMOS[:,1],color='green',yerr=[MeerKAT_COSMOS[:,2],MeerKAT_COSMOS[:,3]],markersize=6,alpha=0.4,label='MIGHTEE COSMOS 1.4GHz (Hale +2022)' ,fmt='^')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_ylim(1e-3, 1e4)
    axs[1].set_xlim(1e-2, 10**2.5)
    axs[1].set_xlabel(r'$S_{1.4 GHz}$ [mJy]', fontsize=16)
    axs[1].set_ylabel(r'$S^{2.5}dN/dS$ [sr$^{-1}$ Jy$^{1.5}$]', fontsize=16)
    axs[1].axvline(5*10e-3,linestyle='--',color='red',alpha=1)
    axs[1].legend(loc='lower right')
    axs[1].set_title('Zwcl2341', fontsize=18)
    axs[1].tick_params(axis='both', which='major', labelsize=12)
    axs[1].tick_params(axis='both', which='minor', labelsize=12)
    plt.tight_layout()
    plt.savefig(f'/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/Source_Counts.png',
                bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.show()
    plt.close()

    # plt.figure(figsize=(10, 7)) 
    # plt.plot(S_M,M_counts,label='Mancuso 2017 SFG + AGN model (Mancuso +2017)',color='black')
    # plt.plot(Semper_M,Semper_counts,label='SEMPER SFG model (SEMPER +2024)',color='green')
    # #plt.plot(TRECS_M,TRECS_count,label='TRECS SFG model (TRECS +2024)',color='orange')
    # plt.plot(S,10**SCs_Bondi(S),label='COSMOS 1.4GHz Bondi fitted (Bondi +2008)')
    # #plt.plot(S_M,M_counts,label='Massardi 2010 SFG + AGN model (Massardi +2010)',color='cyan')
    # #plt.plot(S,10**SCs_CLASS(S*1e-3),label='Super-CLASS Risely fitted (Risely +2016)')
    # #plt.scatter(centres_CLASS*1e3,counts_norm_CLASS,label='SuperCLASS GMRT 325MHz (Risely +2016)',marker='o',color='red')
    # plt.errorbar(MeerKAT_XMM_LSS[:,0]*1e-3,MeerKAT_XMM_LSS[:,1],color="#16A5BE",yerr=[MeerKAT_XMM_LSS[:,2],MeerKAT_XMM_LSS[:,3]],markersize=4,alpha=1,label='MeerKAT XMM LSS' ,fmt='o')
    # plt.errorbar(MeerKAT_COSMOS[:,0]*1e-3,MeerKAT_COSMOS[:,1],color='green',yerr=[MeerKAT_COSMOS[:,2],MeerKAT_COSMOS[:,3]],markersize=4,alpha=1,label='MeerKAT COSMOS' ,fmt='^')
    # plt.errorbar(centres[0]*1e3, counts_norm[0],yerr=counts_err[0],color='#ADD8E6' ,fmt='o',markersize=4,alpha=1,label=f'uncorrected A2631 MeerKAT 1.28GHz (This paper)')
    # plt.errorbar(centres[1]*1e3, counts_norm[1],yerr=counts_err[1],color="#DA6969" ,fmt='o',markersize=4,alpha=1,label=f'uncorrected Zwcl2341 MeerKAT 1.28GHz (This paper)')
    # plt.errorbar(centres[0]*1e3, counts_norm[0]*completeness_corr[0],yerr=counts_err[0],xerr=deltaS[0],color='blue',label=f'A2631 MeerKAT 1.28GHz (This paper)' ,fmt='+',markersize=8,alpha=1)   
    # plt.errorbar(centres[1]*1e3, counts_norm[1]*completeness_corr[1],yerr=counts_err[1],xerr=deltaS[1],color='red',label=f'Zwcl2341 MeerKAT 1.28GHz (This paper)' ,fmt='+',markersize=8,alpha=1) 
    # plt.errorbar(x_lit,y_lit,yerr=[y_err+var,x_err-var], fmt='*',color='purple',label='de Zotti 1.4GHz compilation (de Zotti +2010)',alpha=0.3)
    # plt.errorbar(x_DEEP*1e3,y_DEEP,yerr=[y_err_DEEP,x_err_DEEP], fmt='x',color='orange',label='DEEP2 MeerKAT 1.28GHz  (Mauch +2019)',alpha=0.9)
    # #plt.errorbar(x_LoTSS,y_LoTSS,yerr=[y_LoTSS_err_up, y_LoTSS_err_down], fmt='*',color='orange',label='LoTSS 150MHz  (Mandal +2019)',alpha=0.5)
    # # plt.axvline(5*16e-3,linestyle='--',color='blue',alpha=1)
    # # plt.axvline(5*10e-3,linestyle='--',color='red',alpha=1)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.ylim(10**-3,10**4)
    # plt.xlim(10**-2,10**2.5)
    # plt.xlabel(r'$S_{1.4 GHz}$ [mJy ]',fontsize=18)
    # plt.ylabel(r'$S^{2.5}dN/dS$ [sr$^{-1}$ Jy$^{1.5}$]',fontsize=18)
    # plt.tick_params(axis='both', which='major', labelsize=15, length=3, width=1)  # Increase size of major tick labels
    # plt.tick_params(axis='both', which='minor', labelsize=15, length=3, width=1)
    # plt.legend()
    # plt.title(f'{name}',size=20)
    # plt.tight_layout()
    # plt.savefig(f'/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/Source_Counts_{name}.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    # plt.show()
    # plt.close()

    with open("A2631_source_counts.txt", "w") as file:
        for c, d, n, counts, err, x,y in zip(centres[0],deltaS[0], numbers[0],counts_norm[0], counts_err[0], completeness_corr[0],resolution_corr[0] ):
            file.write(f"{c*1e3} {d*1e3} {n} {counts} {err} {x} {y}  \n")
    with open("Zwcl2341_source_counts.txt", "w") as file:
        for c, d, n, counts, err, x,y in zip(centres[1],deltaS[1], numbers[1],counts_norm[1], counts_err[1], completeness_corr[1],resolution_corr[1] ):
            file.write(f"{c*1e3} {d*1e3} {n} {counts} {err} {x} {y} \n")

if __name__ == "__main__":
    outname='/home/kincaid/Desktop/Saraswati_codes/plots/'
    source_count_literature='de_zotti.txt'
    source_count_literature_DEEP='MeerKAT_deep.txt'
    source_count_literature_SuperCLASS='SuperCLASS_source_counts.txt'
    MeerKAT_COSMOS_lit='MeerKAT_COSMOS.txt'
    MeerKAT_XMM_LSS_lit='MeerKAT_XMM_LSS.txt'
    survey_area=1.65# (0.7*0.73*np.pi) in square degrees
    nbins=25
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
        resolution_bias_file='resolution_interp_func.pkl'
        completeness_file = name+'/'+"visib_correction_cut.txt"
        false_detection_file='false_detection_correction.txt'
        real_catalog =  '/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_eddington_corr_srl.fits'  
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
    corrections(numbers_array,deltaS_array,centres_array,counts_norm_array,counts_err_array,source_count_literature, source_count_literature_DEEP,source_count_literature_SuperCLASS,completeness_corr_array, resolution_corr_array )  #calcaulte the ratio of recovered/injected in each bin
    
    
