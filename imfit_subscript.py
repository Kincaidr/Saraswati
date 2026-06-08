import sys

def find_flux(image, region,index):
    print('image',image)
    print('line',index)
    lines= open(region).readlines()
    try:
        MH_output=imfit(imagename=image,region=lines[index])
        flux=MH_output['results']['component0']['flux']['value'][0]*1e3
        flux_err=MH_output['results']['component0']['flux']['error'][0]*1e3
    except KeyError:
                flux=0
                flux_err=0
    return(flux,flux_err)

if __name__== "__main__":
    image = sys.argv[1]
    region = sys.argv[2]
    index = int(sys.argv[3])
    flux,flux_err=find_flux(image,region,index)


    