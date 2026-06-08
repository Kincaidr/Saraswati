from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle


def position(RA_array,DEC_array,RA_full_array,DEC_full_array):
    cluster_centre1=SkyCoord(str(354.41916666667 ), str(0.27666666666667), frame='icrs',unit=(u.deg,u.deg))
    cluster_centre2=SkyCoord(str(355.91541666667), str(0.33083333333333), frame='icrs',unit=(u.deg,u.deg))

    plt.figure(figsize=(8, 6)) 
    plt.scatter(RA_array[0],DEC_array[0],color='#56B4E9',s=5,alpha=0.4, label='A2631')
    plt.scatter(RA_array[1],DEC_array[1],color='#D55E00',s=5,alpha=0.4,label='ZwCL2341')
    plt.scatter(RA_full_array[0],DEC_full_array[0],color='blue',s=4,alpha=0.2)
    plt.scatter(RA_full_array[1],DEC_full_array[1],color='red',s=4,alpha=0.2)
    plt.plot(cluster_centre1.ra, cluster_centre1.dec, marker='+', color='black', markersize=35, markeredgewidth=3)
    plt.plot(cluster_centre2.ra, cluster_centre2.dec, marker='+', color='black', markersize=35, markeredgewidth=3)
    plt.xlabel('RA [deg]',size=20)
    plt.ylabel('DEC [deg]',size=20)
    x0 = 353.2  # starting RA position for the scale bar
    y0 = min(DEC_array[0]) - 0.4
    plt.plot([x0, x0 + 1], [y0, y0], color='black', linewidth=3)  # 1-degree long bar
    plt.text(x0 + 0.5, y0 + 0.15, '1°', ha='center', va='top', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=15, length=3, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=15, length=3, width=1)
    plt.legend(fontsize=16)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/RA_DEC_positions.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()


if "__main__":
    z=0.27
    names=['A2631','Zwcl2341']
    RA_array=[]
    DEC_array=[]
    RA_full_array=[]
    DEC_full_array=[]
    for name in names:
        path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
        masked_cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_cut_srl.fits')
        full_cat=Table.read('/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_srl.fits')
        RA_array.append(masked_cat['RA'])
        DEC_array.append(masked_cat['DEC'])
        RA_full_array.append(full_cat['RA'])
        DEC_full_array.append(full_cat['DEC'])
    position(RA_array,DEC_array,RA_full_array,DEC_full_array)


