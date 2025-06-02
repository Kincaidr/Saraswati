from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def get_area(fits_image):
    hdu = fits.open(fits_image)[0]
    data = hdu.data
    if data.ndim > 2:
        data = data[0, 0, :, :]
    rmsmin = np.linspace(0, np.nanmax(data), 501)
    areas = np.zeros_like(rmsmin)
    for i, r in enumerate(rmsmin):
        areas[i] = data[data <= r].size * (1.5 / 3600**2)
    full_area = data[~np.isnan(data)].size * (1.5 / 3600**2)
    return(areas,full_area,rmsmin)

def plot(areas,full_areas, rmsmin):
    fig, ax1 = plt.subplots(figsize=(8, 8))
    # Original plot (blue)
    ax1.plot(rmsmin[0] * 1e3, areas[0], color='blue',label='A2631')
    ax1.plot(rmsmin[1] * 1e3, areas[1], color='red',label='Zwcl2341')
    ax1.set_xlabel(r'RMS (mJy/Beam)', size=20)
    ax1.set_ylabel('Area (deg²)', size=20)
    ax1.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=14)
    ax1.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

    ax1.set_xlim(0,1)
    ax2 = ax1.twinx()
    ax2.set_ylim((ax1.get_ylim()[0] /full_areas[0])*100, ax1.get_ylim()[1]/full_areas[0] * 100)
    ax2.set_ylabel('Area Percent (%)', size=20)
    ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
    ax2.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)
    ax1.legend(loc='upper left', fontsize=14) 
    plt.tight_layout()
    plt.savefig('plots/Visibility_function.png',bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()


if "__main__":
    names=['A2631','Zwcl2341']
    Areas=[]
    rmsmins=[]
    full_areas=[]
    for name in names:
        fits_image= '/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'+name+'_full_rms_map.fits'
        areas,full_area,rmsmin=get_area(fits_image)
        Areas.append(areas)
        rmsmins.append(rmsmin)
        full_areas.append(full_area)
    plot(Areas,full_areas,rmsmins)