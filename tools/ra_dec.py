#


from astropy.io import fits
from astropy.table import Table
import sys
catalog=sys.argv[1]
cat=Table.read(catalog)

RA, DEC=cat['RA'], cat['DEC']
new_table = Table([RA, DEC], names=('RA', 'DEC'))

# Save the new table to a file (e.g., FITS format)
new_table.write('ra_dec_table.fits', format='fits', overwrite=True)
print('new table written')
