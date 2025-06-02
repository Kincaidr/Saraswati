import numpy as np
import sys

image1=sys.argv[1]
frequency1=int(sys.argv[2])
image2=sys.argv[3]
frequency2=int(sys.argv[4])
region='A2631_VLA_circle_regions.crtf'
print('region',region)
print('freq1 is', frequency1)
print('freq2 is', frequency2)
lines = open(region).readlines()
print(lines[1])
flux_values1 = []
flux_values2 = []
flux_values1_err = []
flux_values2_err = []


 # Create a new list for each image
def linear_function(x, m, b):
    return m * x + b

for j in range(2, len(lines)):
    image=image1
    MH_output=imfit(imagename=image,region=lines[j])
    flux=MH_output['results']['component0']['flux']['value'][0]*1e3
    flux_err=MH_output['results']['component0']['flux']['error'][0]*1e3
    flux_values1.append(flux) #mJy
    flux_values1_err.append(flux_err)

for j in range(2, len(lines)):
    image=image2
    MH_output=imfit(imagename=image,region=lines[j])
    flux=MH_output['results']['component0']['flux']['value'][0]*1e3
    flux_err=MH_output['results']['component0']['flux']['error'][0]*1e3
    flux_values2.append(flux) #mJy
    flux_values2_err.append(flux_err)

print('Flux values min',np.min(flux_values1) )
print('Flux values max',np.max(flux_values2))
flux_values1 = np.array(flux_values1)
flux_values2 = np.array(flux_values2)
flux_values2_err = np.array(flux_values2_err)
flux_values1_err = np.array(flux_values1_err)
#frequency1=1.594e+09
#frequency2=1.283e+09

if frequency1 > frequency2:
    flux_values2_corr= flux_values2* (frequency1 / frequency2)**-0.7
    flux_values2_corr_err= flux_values2_err * (frequency1 / frequency2)**-0.7

    flux_values1_corr = flux_values1
    flux_values1_corr_err= flux_values1_err 
    print(flux_values1_corr,flux_values2_corr)
else:  
    flux_values1_corr= flux_values1* (frequency2 / frequency1)**-0.7
    flux_values1_corr_err= flux_values1_err * (frequency1 / frequency2)**-0.7
    
    flux_values2_corr = flux_values2
    flux_values2_corr_err= flux_values2_err 
    print(flux_values1_corr,flux_values2_corr)

path="/home/kincaid/Desktop/Saraswati_codes/A2631/flux_scale/wsclean_VLA/"
with open(path+"Flux_comparison.txt", "w") as file:
    for flux1, flux1_err, flux2, flux2_err in zip(flux_values1_corr,flux_values1_corr_err,  flux_values2_corr, flux_values2_corr_err):
        file.write(f"{flux1} {flux1_err} {flux2}  {flux2_err} \n")


