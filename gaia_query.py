
from astropy.table import Table, vstack
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u

breakpoint()
slice_size=5000

#cluster_centre=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
cluster_centre=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))
radius=1*u.deg
offset=0
#sql_query=f""" SELECT b.ra, b.dec, b.photo_g_mean_flux  FROM gaiadr3.gaia_source as b WHERE Q3C_RADIAL_QUERY(b.ra,b.dec,{cluster_centre.ra.value}, {cluster_centre.dec.value},2)"""

#sql_query=f""" SELECT b.ra, b.dec  FROM gaia_dr2.gaia_source as b WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{cluster_centre.ra.value}, {cluster_centre.dec.value},1))"""

result = Vizier(columns=["RA_ICRS", "DE_ICRS"]).query_region(
    cluster_centre,
    radius=radius,
    catalog="I/355/gaiadr3"
)

gaia_table = result[0]
print(len(gaia_table))
gaia_table.write("gaia_vizier_result.fits", format="fits", overwrite=True)



# total_rows_str = qc.query(sql=sql_query, scalar=True)
# result = total_rows_str.split('\n')[1:-1]
# total_rows=len(result)
# print('total rows',total_rows)
# combined_result = []

# # Iterate through the data in smaller slices
# for offset in range(0, total_rows, slice_size):
#     t1 = time.time()
#     # Use qc.query() to execute the query and get the result
#     result = qc.query(sql=sql_query)
#     result1 = result.split('\n')[1:-1]
#     t2 = time.time()
#     print('Time',t2-t1, len(result1), slice_size, offset)
  
#     # If the result is empty, it means we have processed all rows, so break the loop
#     if not result:
#         break
#     else:
#         keys = np.array(result.split('\n')[0].split(','))
#         t2 = time.time()
#         print('Query finished in ', t2-t1, 'seconds')
#         data=[]
#         for row_str in result1:
#             values = row_str.split(',')
#             data.append(dict(zip(keys, values)))

#         table = Table(rows=data)
#         output_filename = 'query_result_'+str(offset)+'.fits'
#         table.write(output_filename, format='fits', overwrite=True)
#         print('table '+ output_filename + ' has been written ')


