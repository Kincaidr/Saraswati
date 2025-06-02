import os
import shutil

number1=(sys.argv[1])
number2=(sys.argv[2])

importfits(fitsimage=number1,imagename=number1+'.img',overwrite=True)
importfits(fitsimage=number2,imagename=number2+'.img',overwrite=True)

mode='get'

cdet1_1=imhead(imagename=number1+'.img',mode=mode,hdkey='CDELT1')
cdet1_2=imhead(imagename=number1+'.img',mode=mode,hdkey='CDELT2')
cdet2_1=imhead(imagename=number2+'.img',mode=mode,hdkey='CDELT1')
cdet2_2=imhead(imagename=number2+'.img',mode=mode,hdkey='CDELT2')

cellsize1_1=cdet1_1['value']*57.2958*3600
cellsize1_2=cdet1_2['value']*57.2958*3600

cellsize2_1=cdet2_1['value']*57.2958*3600
cellsize2_2=cdet2_2['value']*57.2958*3600

if abs(cellsize1_1)<abs(cellsize2_1) and cellsize1_2<cellsize2_2:
 temp=str(number1+'.img');pbc=(number2+'.img')
else:
 temp=str(number2+'.img');pbc=(number1+'.img') 
 
out3=pbc+'.rgd'

imregrid(imagename=pbc,template=temp,asvelocity=False,output=out3)
exportfits(imagename=out3,fitsimage=out3+'.fits', overwrite=True)

bmj1_1=imhead(imagename=pbc,mode=mode,hdkey='BMAJ')
bmm1_2=imhead(imagename=pbc,mode=mode,hdkey='BMIN')

beamsize1_1=bmj1_1['value']
beamsize1_2=bmm1_2['value']

bmj2_1=imhead(imagename=temp,mode=mode,hdkey='BMAJ')
bmm2_2=imhead(imagename=temp,mode=mode,hdkey='BMIN')

beamsize2_1=bmj2_1['value']
beamsize2_2=bmm2_2['value']

if beamsize1_1 > beamsize2_1 and beamsize1_2>beamsize2_2:
 x=round(beamsize1_1);y=round(beamsize1_2);fitsimage=temp
else:  
 x=round(beamsize2_1);y=round(beamsize2_2);fitsimage=out3

z=0

major=str(x)+'arcsec'
minor=str(y)+'arcsec'
pa=str(z)+'deg'


outfile=fitsimage+'.conv'
image2=fitsimage+'.conv.fits'

imsmooth(imagename=fitsimage,targetres=True,major=major,minor=minor,pa=pa,outfile=outfile,overwrite=True)
exportfits(imagename=outfile,fitsimage=image2, overwrite=True)

#shutil.rmtree(pbc)
#shutil.rmtree(temp)
#shutil.rmtree(out3)
#shutil.rmtree(outfile)

