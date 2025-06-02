
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def create_output_file_v1(output_filename, visiblity_correction_file,inj_sources, no_simulations, plots, name ):

    data = np.loadtxt(output_filename )
    rec_counts = np.sum(data[1:11, :], axis=0)
    data_corr= np.loadtxt(visiblity_correction_file)
    area=data_corr[:,0]
    bin_centres=data_corr[:,1]
    bin_extrap=data[0,:]
    rec_counts=np.array(rec_counts)
    rec_counts_corr=np.array(rec_counts)*area

    for i in range(len(rec_counts)):
        if rec_counts[i] > inj_sources*no_simulations:
            rec_counts[i] = inj_sources*no_simulations

    print("Bin Centres:", bin_centres)
    print("Summed Recovered Counts:", rec_counts)

    inj_sources_tot=np.full(nbins,inj_sources*no_simulations)
    ratio = rec_counts /inj_sources_tot
    error=ratio*np.sqrt((np.sqrt(inj_sources_tot)/inj_sources_tot)**2+(np.sqrt(rec_counts)/rec_counts)**2)
    ratio_corr = rec_counts_corr /inj_sources_tot
    error_corr=ratio*np.sqrt((np.sqrt(inj_sources_tot)/inj_sources_tot)**2+(np.sqrt(rec_counts_corr)/rec_counts_corr)**2)

    output_file = name+"_output_table.txt"
    with open(output_file, "w") as file:
        for r, value, value_err, flux, b in zip(ratio, ratio_corr, error_corr, bin_centres, bin_extrap):
            file.write(f"{r} {value} {value_err} {flux} {b}\n")

    print(f"Ratio array has been written to {output_file}")

def create_output_file(output_filename,inj_sources, no_simulations, name ):

    data = np.loadtxt(output_filename )
    rec_counts = np.sum(data[1:11, :], axis=0)
    bin_centres = data[0, :]

    for i in range(len(rec_counts)):
        if rec_counts[i] > inj_sources*no_simulations:
            rec_counts[i] = inj_sources*no_simulations

    print("Bin Centres:", bin_centres)
    print("Summed Recovered Counts:", rec_counts)

    inj_sources_tot=np.full(nbins,inj_sources*no_simulations)
    ratio = rec_counts /inj_sources_tot
    
    error=np.sqrt((np.sqrt(inj_sources_tot)/inj_sources_tot)**2+(np.sqrt(rec_counts)/rec_counts)**2)
    output_file = name+"_output_table.txt"
    with open(output_file, "w") as file:
        for r, rat, err in zip(bin_centres,ratio, error ):
            file.write(f"{r} {rat} {err}\n")

    print(f"Ratio array has been written to {output_file}")

    return(output_file)

def interpolation(visiblity_correction_file,output_filename):                                    

    data=np.loadtxt(output_filename)
    extrap_data=np.loadtxt(visiblity_correction_file)
    flux_extrap=extrap_data[:,1]
    ratio_extrap=extrap_data[:,0]
    bin_centres=data[:,0]

    f=interp1d(flux_extrap, ratio_extrap, bounds_error=False, fill_value="extrapolate")
    correction_factor=f(bin_centres)
    breakpoint()
    return(bin_centres,correction_factor)

def plot(output_filename,bin_centres, corr):

    data=np.loadtxt(output_filename)
    comp=data[:,1]
    comp_err=data[:,2]
    x=np.array(bin_centres)
    completeness=np.array(comp)
    completeness_corr=np.array(corr)
    completeness_err=np.array(comp_err)
    breakpoint()
    #err_corr=np.array(error_co)
    plt.errorbar(x,completeness,yerr=completeness_err, label='Fraction')
    plt.errorbar(x,completeness_corr*completeness,yerr=completeness_err, label='Visiblity corrected fraction')
    plt.xlabel(r'Log Flux $S_T$ [mJy]',size=18)
    plt.ylabel('Detected Fraction',size=18)
    plt.tick_params(axis='both', which='major', labelsize=15, length=5, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=15, length=5, width=1)
    #plt.axvline(x=sigma, color='red', linestyle='--', label=r'detection    threshold $4 \sigma$')
    plt.axhline(y=1,color='green')
    plt.xscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plots+name+'_ratio_plot.png')
    plt.show()
    plt.close()


if "__main__":
    name='A2631'
    inj_sources=500
    no_simulations=20
    nbins=30
    completeness_filename = name+"_recovered_counts_table.txt"
    visiblity_correction_file=name+'_corrected_area.txt'
    plots='/home/kincaid/Desktop/Saraswati_codes/'+name+'/plots/'
    completeness_filename=create_output_file(completeness_filename,inj_sources, no_simulations,name )
    bin_centre, corr=interpolation(visiblity_correction_file,completeness_filename)    
    plot(completeness_filename,bin_centre,corr)