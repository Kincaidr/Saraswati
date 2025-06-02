import sys
import numpy as np

region='20_regions.crtf'
image='simulated_image_test_n0.fits'

print('region is ', region)

print('image is ', image)

lines = open(region).readlines()
flux=np.zeros(len(lines))

for i in range(len(lines)):
    #f = open("fluxes.txt", "aw")

    fit=imfit(imagename=image,region=lines[i])
    comp=fit['results']['component0']
    #breakpoint()
    flux[i]=comp['flux']['value'][0]

        #f.write(str(flux)+'\n')
        







   
 
 
 
 


