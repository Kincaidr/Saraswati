from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
#cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
radius=1*u.deg
offset=0
#sql_query=f""" SELECT b.ra, b.dec, b.photo_g_mean_flux  FROM gaiadr3.gaia_source as b WHERE Q3C_RADIAL_QUERY(b.ra,b.dec,{cluster_centre.ra.value}, {cluster_centre.dec.value},2)"""

#sql_query=f""" SELECT b.ra, b.dec  FROM gaia_dr2.gaia_source as b WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{cluster_centre.ra.value}, {cluster_centre.dec.value},1))"""

vizier = Vizier()
vizier.ROW_LIMIT = -1
result = vizier.query_region(cluster_centre,radius=radius, catalog="VIII/59/first")
table = result[0]
RA,DEC=table['RAJ2000'],table['DEJ2000']    
breakpoint()
coord = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg), frame='icrs')
ra_deg = coord.ra.deg
dec_deg = coord.dec.deg
table['RA'] = ra_deg
table['DEC'] = dec_deg
print(len(table))
table.write("A2631_FIRST_vizier_result.fits", format="fits", overwrite=True)

