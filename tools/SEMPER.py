import pickle
import numpy as np
from scipy.interpolate import interp1d

def SEMPER_SFG_AGN_counts():
    file1='tools/AGN_model.txt'
    file2='tools/number_counts_SFGs_SEMPER.pkl'
    AGN_model=np.loadtxt(file1)
    S_AGN=np.log10((10**AGN_model[:,0])*1e3) # Convert to mJy
    FSRQ=AGN_model[:,1]
    BL_Lac=AGN_model[:,2]
    steep_spectrum=AGN_model[:,3]

    with open(file2, "rb") as file2:
        data = pickle.load(file2)

    new_grid = np.linspace(-3.9, 5.8, 10000) 
    FSRQ_func = interp1d(S_AGN, FSRQ, bounds_error=False, fill_value=0)
    BL_Lac_func = interp1d(S_AGN, BL_Lac, bounds_error=False, fill_value=0)
    steep_spectrum_func = interp1d(S_AGN, steep_spectrum, bounds_error=False, fill_value=0)
    FSRQ_interp = FSRQ_func(new_grid)
    BL_Lac_interp = BL_Lac_func(new_grid)
    steep_spectrum_interp = steep_spectrum_func(new_grid)

    AGN_counts=  FSRQ_interp + BL_Lac_interp + steep_spectrum_interp
    counts_semper=data(new_grid)
    total =AGN_counts+10**counts_semper
    return 10**new_grid, total

