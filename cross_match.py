import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.table import Table

def cross_match_nearest_neighbor(tab1_path, tab2_path):
    # Read in the catalogs
    table1 = Table.read(tab1_path)
    table2 = Table.read(tab2_path)
    print(f'Length of table1: {len(table1)}')
    print(f'Length of table2: {len(table2)}')
    coords1 = SkyCoord(ra=table1['RA'], dec=table1['DEC'])
    coords2 = SkyCoord(ra=table2['ra']*u.deg, dec=table2['dec']*u.deg)
    idx, d2d, _ = match_coordinates_sky(coords1, coords2)
    d2d_arcsec = d2d.arcsec
    mask=d2d_arcsec < 8
    d2d_arcsec = d2d_arcsec[mask]
    return d2d_arcsec  # Return if you want to analyze further

# Usage
if __name__ == "__main__":
    clusters = ['A2631', 'Zwcl2341']
    path = f'/home/kincaid/Desktop/Saraswati_codes/catalogs/'
    tab1 = path + 'MeerKAT_combined.fits'
    tab2 = path + 'A2631_HSC_dnnzi.fits'
    colors=['blue', 'red']
    fig, ax = plt.subplots(figsize=(8, 6))
    for i,cluster in enumerate(clusters):
        tab2 = path + f'{cluster}_HSC_dnnzi.fits'
        print(f'Cross-matching {tab1} with {tab2}')
        d2d_arcsec = cross_match_nearest_neighbor(tab1, tab2)

        # Plot histogram of separations for each cluster on the same axes
        ax.hist(d2d_arcsec, bins=50, histtype='step',color=colors[i], label=cluster)

    ax.set_xlabel('Separation (arcsec)',fontsize=17)
    ax.set_ylabel('Number',fontsize=17)
    ax.axvline(3, color='black', linestyle='--', linewidth=3)
    ax.grid(True)
    ax.legend(fontsize=17)
    ax.tick_params(axis='both', which='major', labelsize=16, length=5, width=1)  # Increase size of major tick labels
    ax.tick_params(axis='both', which='minor', labelsize=16, length=5, width=1)
    plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/HSC_crossmatch.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.show()
