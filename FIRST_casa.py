import glob
import re

def write_file(flux_values,err_values):     
    with open("First_flux.txt", "w") as file:
        for flux, err, in zip(flux_values,err_values):
            file.write(f"{flux} {err}\n")

def find_flux(image,line):
    print('image',image)
    print('line',line)

    try:
        MH_output=imfit(imagename=image,region=line)
        flux=MH_output['results']['component0']['flux']['value'][0]*1e3
        flux_err=MH_output['results']['component0']['flux']['error'][0]*1e3
    except KeyError:
                flux=0
                flux_err=0
    return(flux,flux_err)

if __name__== "__main__":
    region='A2631/A2631_VLA_circle_regions.crtf'
    lines= open(region).readlines()
    path='/home/kincaid/Desktop/Saraswati_codes/A2631/images/cutouts/'
    image_cutouts=glob.glob(path+'first_cutout_*_convolved.fits')
    image_cutouts.sort(key=lambda x: int(re.findall(r'first_cutout_(\d+)_convolved\.fits', x)[0]))
    flux_values= []
    flux_values_err = []
    for image,line in zip(image_cutouts,lines):
        flux_value,flux_err=find_flux(image,line)
        flux_values.append(flux_value)
        flux_values_err.append(flux_err)
    write_file(flux_values,flux_values_err)