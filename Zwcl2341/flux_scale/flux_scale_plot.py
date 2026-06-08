
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

def frequency_scaling(flux, frequency1, frequency2, spectral_index):
    return flux * (frequency2 / frequency1) ** spectral_index

def flux_plot(data):
    flux_values1 = data[:,0]
    flux_values1_err = data[:,1]
    flux_values2 = data[:,2]
    flux_values2 = frequency_scaling(flux_values2, 1.283, 3, -0.7)
    flux_values2_err = data[:,3]
    flux_values2_err = frequency_scaling(flux_values2_err, 1.283, 3, -0.7)

    #slope, intercept,r_value, p_value, std_err = stats.linregress(np.log10(flux_values1), np.log10(flux_values2))
    xx = np.linspace(0,50,1024)
    fig = plt.figure(figsize=(13, 10))

    plt.plot(xx,xx, "k--")
    #plt.plot(xx,linear_function(xx, slope, intercept),'-',linewidth=2.0, label=r'$\alpha$'+' = ' '%0.2f' % (r_value),color='orange')
    #plt.scatter(flux_values1_corr ,flux_values2_corr,label=r'Differences % = ' + '%0.2f' % np.abs(err.mean()))
    plt.errorbar(flux_values1, flux_values2,xerr=flux_values1_err,yerr=flux_values2_err,fmt='+', markersize=15)
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('MeerKAT Flux Density [mJy]',size=18)
    plt.xlabel('VLA Flux Density [mJy]',size=18)
    plt.tick_params(axis='both', which='major', labelsize=18, length=5, width=1)  # Increase size of major tick labels
    plt.tick_params(axis='both', which='minor', labelsize=18, length=5, width=1)
    plt.savefig('flux_scale.png')
    plt.show()


if __name__ == "__main__":

    plots="/home/kincaid/Desktop/Saraswati_codes/A2631/plots/"
    data=np.loadtxt('Flux_comparison.txt')
    flux_plot(data)
