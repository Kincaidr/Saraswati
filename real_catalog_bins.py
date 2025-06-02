import numpy as np
from astropy.table import Table



real_catalog = '/home/kincaid/Desktop/Saraswati_codes/A2631_new/A2631_new_srl.fits'
rec_cat = Table.read(real_catalog)

flux=rec_cat['Total_flux']
nbins=30
pers = np.linspace(0, 100, nbins+1)
Range_x = np.percentile(flux, pers)
breakpoint()
centres = (Range_x[0:-1] + Range_x[1:]) / 2.0

centres=np.array(centres)

output_file = "real_catalog_bins.txt"
with open(output_file, "w") as file:
    for value in centres:
        file.write(f"{value}\n")

print(f"real catalog bins has been written to {output_file}")