
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def plot(real_cat, inverted_cat, plots):
    cat = Table.read(real_cat)
    inv_cat=Table.read(inverted_cat)
    flux=cat['Total_flux']*1e3
    rms=inv_cat['Isl_rms']*1e3
    flux_inv=inv_cat['Total_flux']*1e3
    SNR=flux_inv/rms
    print('Total number of false detections',len(flux_inv))
    print('Total number of false detections SNR > 10 ',len(SNR[SNR>10]))
    print('Fraction of spurious sources',len(flux_inv)/len(flux)*100)
    nbins=27
    Range_x = 10**np.linspace(start=np.log10(flux.min()), stop=np.log10(50), num=nbins+1)
    centres = (Range_x[0:-1] + Range_x[1:]) / 2.0
    hist, x = np.histogram(flux, bins=Range_x)
    hist_inv, x = np.histogram(flux_inv, bins=Range_x)
    ratio=hist_inv/hist
    correction=1-ratio
    with open("false_detection_correction.txt", "w") as file:
        for a,b in zip(centres,correction):
            file.write(f"{a} {b}  \n")

    # plt.figure(figsize=(8, 6))
    # plt.plot(centres, correction, color='blue', label='False detection correction')
    # plt.xlabel(r"Flux [mJy]", size=20)
    # plt.ylabel(r"Counts", size=20)
    # plt.tick_params(axis='both', which='both', direction='in', length=6 , labelsize=14)  # 'in' means they point inward
    # plt.tick_params(which='minor', length=3)
    # plt.xscale('log')
    # plt.savefig(plots+"False_detections.png")
    # plt.show()


if "__main__":
    names=['A2631','Zwcl2341']
    for name in names:
        real_cat = '/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_cut_srl.fits'
        inverted_cat='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'+name+'_cut_inv_srl.fits'
        plots='/home/kincaid/Desktop/Saraswati_codes/'+name+ '/plots/'
        plot(real_cat, inverted_cat, plots)
        

