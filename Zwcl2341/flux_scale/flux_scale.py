import numpy as np
import sys

image1=sys.argv[1]
image2=sys.argv[2]
region='zwcl_regions_2.crtf'
print('region',region)
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
    print('number', j)
    print('region is', lines[j])
    MH_output=imfit(imagename=image,region=lines[j]) 
    flux=MH_output['results']['component0']['flux']['value'][0]*1e3
    flux_err=MH_output['results']['component0']['flux']['error'][0]*1e3
    flux_values1.append(flux) #mJy
    flux_values1_err.append(flux_err)

for j in range(2, len(lines)):
    image=image2
    print('number', j)
    print('region is', lines[j])
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

path="/home/kincaid/Desktop/Saraswati_codes/Zwcl2341/flux_scale/"
with open(path+"Flux_comparison.txt", "w") as file:
    for flux1, flux1_err, flux2, flux2_err in zip(flux_values1,flux_values1_err,  flux_values2, flux_values2_err):
        file.write(f"{flux1} {flux1_err} {flux2}  {flux2_err} \n")


