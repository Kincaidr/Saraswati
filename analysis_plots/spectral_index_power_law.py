import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import json

def find_flux(image,freq):
    flux_values = []
    freq_values = []
    for j in range(2, len(lines)):
        print('Source',j)
        MH_output=imfit(imagename=image,region=lines[j])
        flux=MH_output['results']['component0']['flux']['value'][0]
        flux_values.append(flux) #mJy
        freq_values.append(freq)
    return(flux_values,freq_values)
    
def plot(json_filename, images):
    with open(json_filename, "r") as json_file:
        dict = json.load(json_file)
        
    SI=[]
    for i in range(len(lines)-2):
        fluxes=[]
        freqs=[]
        for image in images:
            flux=dict[image]['flux'][i]
            freq=dict[image]['freq'][i]
            fluxes.append(flux)
            freqs.append(freq)
     
        mask=np.array(fluxes) >0
        fluxes=np.array(fluxes)[mask]
        freqs=np.array(freqs)[mask]
        Result = linregress(np.log10(freqs), np.log10(fluxes))
        SI.append(Result[0])

    plt.figure(figsize=(12, 8))
    path='/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'
    #plt.subplot(1, 2, 1)
    plt.hist(SI, bins=30, color='blue', alpha=0.7)
    plt.xlabel("Spectral Index (α)",size=20)
    plt.ylabel("Number of Sources",size=20)
    #plt.title("Spectral Index Distribution")
    plt.axvline(np.median(SI), color='red', linestyle='dashed', label=f'Median: {np.median(SI):.2f}')
    plt.legend(fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15, length=5, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=15, length=5, width=1)
    plt.savefig(path+'Spectral_Index_Distribution.png')
    plt.show()

if "__main__":

    json_filename = "flux_freq_data.json"
    region='/home/kincaid/Desktop/Saraswati_codes/A2631/flux_scale/A2631_spectral_regions.crtf'
    path='/home/kincaid/Desktop/Saraswati_codes/A2631/spectral/'
    images=["mypipelinerun_ABELL2631_4-0000-image.fits","mypipelinerun_ABELL2631_4-0001-image.fits","mypipelinerun_ABELL2631_4-0002-image.fits","mypipelinerun_ABELL2631_4-0003-image.fits"]#,"mypipelinerun_ABELL2631_4-0004-image.fits"]
    freqs=[941,1111,1283]
    lines = open(region).readlines()

    # flux_freq_dict={}
    # for freq,image in zip(freqs,images):
    #     flux_values, freq_values=find_flux(path+image, freq)
    #     flux_freq_dict[image] = {"flux": flux_values, "freq": freq_values} 
    #     with open(json_filename, "w") as json_file:
    #         json.dump(flux_freq_dict, json_file, indent=4)

    plot(json_filename, images)


# Save dictionary to a JSON file
