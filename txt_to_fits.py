import numpy as np
from astropy.io import fits

cat='Zwcl2341_SIMBAD.txt'
ra_list = []
dec_list = []
otype_list = []

with open(cat, 'r') as f:  # <-- replace 'yourfile.txt' with your actual filename
    for line in f:
        parts = line.strip().split()
        if len(parts) == 3:
            ra = float(parts[0])
            dec = float(parts[1])
            otype = parts[2].strip('"') 
            ra_list.append(ra)
            dec_list.append(dec)
            otype_list.append(otype)

# Create FITS columns
col1 = fits.Column(name='RA', format='D', array=np.array(ra_list))
col2 = fits.Column(name='DEC', format='D', array=np.array(dec_list))
col3 = fits.Column(name='OTYPE', format='20A', array=np.array(otype_list))


cols = fits.ColDefs([col1, col2, col3])
hdu = fits.BinTableHDU.from_columns(cols)

# Save to a FITS file
hdu.writeto('Zwcl2341_SIMBAD.fits', overwrite=True)


breakpoint()
