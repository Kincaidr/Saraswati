from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

# Return all matching rows
Vizier.ROW_LIMIT = -1
name='A2631'
# Define coordinates and 60 arcsecond search radius
coord = SkyCoord(ra=354.4191666, dec=0.27666666, unit=(u.deg, u.deg), frame='icrs')
radius = 60 * u.arcsec

# Example: HSC PDR2 catalog in Vizier
# This catalog has photometric redshifts and multi-band photometry
catalog_id = "II/342/hsc2"

# Query the catalog
result = Vizier.query_region(coord, radius=radius, catalog=catalog_id)

breakpoint()
# Access the data table
data = result[0]

data_table = result[0]

# Save as FITS file
output_filename = f"{name}_hsc_photoz_query.fits"
data_table.write(output_filename, format="fits", overwrite=True)

print(f"Saved {len(data_table)} rows to '{output_filename}'")

# Show some rows
print(data[:5])