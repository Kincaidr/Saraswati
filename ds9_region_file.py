from astropy.table import Table

name='A2631'  # Change this to the desired cluster name
path=f"/home/kincaid/Desktop/Saraswati_codes/{name}/catalogs/"

catalog_file = path+name+"_srl.fits"  # e.g., a source catalog
table = Table.read(catalog_file)

ra_col = "RA"     # Change if your RA column has a different name
dec_col = "DEC"   # Change if your Dec column has a different name
shape = "circle"  # Could be 'circle', 'box', etc.
radius_arcsec = 10  # Radius in arcseconds for circle regions

region_lines = [
    "# Region file format: DS9 version 4.1",
    "global color=green font=\"helvetica 10\"",
    "fk5"]

for row in table:
    ra = row[ra_col]
    dec = row[dec_col]
    region_lines.append(f"{shape}({ra},{dec},{radius_arcsec}\")")

with open("catalog_ds9_regions.reg", "w") as f:
    f.write("\n".join(region_lines))

print(f"DS9 region file written:")