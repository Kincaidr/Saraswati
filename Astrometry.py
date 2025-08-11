        
from astropy import coordinates as coords
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
import numpy as np

def astrometry_VLASS(MeerKAT_NVSS_cat):
        MeerKAT_NVSS_cat_fits = fits.open(MeerKAT_NVSS_cat)
        MeerKAT_NVSS_cat=MeerKAT_NVSS_cat_fits[1].data

        MeerKAT_coord=SkyCoord(MeerKAT_NVSS_cat['RA']*u.deg,MeerKAT_NVSS_cat['DEC']*u.deg, frame='fk5')
        NVSS_coord=SkyCoord(MeerKAT_NVSS_cat['RAJ2000']*u.deg,MeerKAT_NVSS_cat['DEJ2000']*u.deg, frame='fk5')

        #import IPython;IPython.embed()
        
        delta_RA= (MeerKAT_coord.ra.value - NVSS_coord.ra.value)*3600
        delta_DEC= (MeerKAT_coord.dec.value - NVSS_coord.dec.value)*3600
        print('Total number of sources before compact criterion',len(delta_RA))
        sel21 = np.logical_and(MeerKAT_NVSS_cat["Total_flux"] / MeerKAT_NVSS_cat["Peak_flux"] < 1.2,
                       MeerKAT_NVSS_cat["Total_flux"] / MeerKAT_NVSS_cat["Peak_flux"] > 0.8)
        print('Total number of soruces',len(delta_RA[sel21]))
        breakpoint()
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.scatter(delta_DEC[sel21],delta_RA[sel21], c="k", marker="o", alpha=0.3)
       # plt.scatter(delta_DEC[~sel21],delta_RA[~sel21], c="b", marker="o", alpha=0.2)
        plt.xlim(-4,4)
        plt.ylim(-4,4)
        plt.axvline(np.mean(delta_DEC),linestyle='dashed',color='black')
        plt.axhline(np.mean(delta_RA),linestyle='dashed',color='black')

        plt.axvline(np.mean(delta_DEC), color="k", linestyle="--", alpha=0.3,
            label="MeerKAT offset to VLASS $\Delta \\alpha , \Delta\delta$={0:.2f}\",{1:.2f}\"".format(np.mean(delta_RA),np.mean(delta_DEC)))

        plt.axhline(np.mean(delta_RA), color="k", linestyle="--", alpha=0.3)
        circle=plt.Circle((np.mean(delta_DEC), np.mean(delta_RA)), np.std(delta_RA), edgecolor= 'red',facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed')
        plt.gca().add_patch(circle)
        ax.add_patch(circle)
        circle2=plt.Circle((np.mean(delta_DEC), np.mean(delta_RA)), 3*np.std(delta_RA), edgecolor= 'red',facecolor='None', linewidth=2, alpha=1 ,ls = 'dashed')
        plt.gca().add_patch(circle2)
        ax.add_patch(circle2)
        print('mean in RA is',np.mean(delta_RA) )
        print('mean in DEC is',np.mean(delta_DEC) )
        plt.xlabel(r'$\delta_{MeerAKT} - \delta_{VLASS}$ (arcsec)',fontsize=20)
        plt.ylabel(r'$\alpha_{MeerAKT} - \alpha_{VLASS}$ (arcsec)',fontsize=20)
        plt.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=22)
        plt.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=22)
        plt.legend(fontsize=17)
        plt.tight_layout()
        plt.savefig('/home/kincaid/Desktop/MOSS2_paper/paper/Figures/plots/VLASS_astrometry.png', bbox_inches='tight', pad_inches=0.1,dpi=300)
        plt.show()

if "__main__":

    MeerKAT_NVSS_cat='/home/kincaid/Desktop/Saraswati_codes/VLASS_MeerKAT_radio_join.fits'
    output_path='/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'
    astrometry_VLASS(MeerKAT_NVSS_cat)