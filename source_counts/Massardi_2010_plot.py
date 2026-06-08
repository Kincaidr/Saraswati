import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def Massardi_counts():
    SFR_model=np.loadtxt('tools/SFG_model.txt')
    AGN_model=np.loadtxt('tools/AGN_model.txt')
    S_AGN=(10**AGN_model[:,0])*1e3 # Convert to mJy
    FSRQ=AGN_model[:,1]
    BL_Lac=AGN_model[:,2]
    steep_spectrum=AGN_model[:,3]
    S_SFR=(10**SFR_model[:,0])*1e3 # Convert to mJy
    Spirals=10**SFR_model[:,1]
    Starburst=10**SFR_model[:,2]
    FSRQ_func = interp1d(S_AGN, FSRQ, bounds_error=False, fill_value=0)
    BL_Lac_func = interp1d(S_AGN, BL_Lac, bounds_error=False, fill_value=0)
    steep_spectrum_func = interp1d(S_AGN, steep_spectrum, bounds_error=False, fill_value=0)
    Spirals_func = interp1d(S_SFR, Spirals, bounds_error=False, fill_value=0)
    Starburst_func = interp1d(S_SFR, Starburst, bounds_error=False, fill_value=0)
    # Define the new flux density grid for plotting and interpolation
    new_grid = np.linspace(1e-3, 200, 10000000)  # 0.01 mJy to 100 mJy
    FSRQ_interp = FSRQ_func(new_grid)
    BL_Lac_interp = BL_Lac_func(new_grid)
    steep_spectrum_interp = steep_spectrum_func(new_grid)
    Spirals_interp = Spirals_func(new_grid)
    Starburst_interp = Starburst_func(new_grid)  
    total_N=FSRQ_interp + BL_Lac_interp + steep_spectrum_interp + Spirals_interp + Starburst_interp
    return(new_grid,total_N)
  

def plot_source_counts(new_grid,total_N):
    # plt.plot(S_AGN, FSRQ, label='FSRQ', color='blue')
    # plt.plot(S_AGN, BL_Lac, label='BL Lac', color='red')
    # plt.plot(S_AGN, steep_spectrum, label='Steep Spectrum', color='green')
    # plt.plot(S_SFR, Spirals, label='Spirals', color='orange')
    # plt.plot(S_SFR, Starburst, label='Starburst', color='purple')
    plt.plot(new_grid, total_N, label='All', color='cyan')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Flux Density (mJy)')
    plt.ylabel('Number of Sources')
    plt.title('Source Counts for AGN and SFR Models')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    new_grid, total_N = Massardi_counts()
    plot_source_counts(new_grid, total_N)
