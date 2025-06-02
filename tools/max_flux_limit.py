import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import root_scalar
from scipy.optimize import brentq

# Assuming the following are defined:
# - SCs_Bondi(Svals) : Source count curve function
# - survey_area : Survey area in square degrees
# - min_flux : Lower flux limit
# - target_num_sources : Total number of sources you want in the image
#extra
def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def SCs_Mandal(S,  a0=1.655, a1=-0.1150, a2=0.2272, a3=0.51788, a4=-0.449661, a5=0.160265, a6=-0.028541, a7=0.002041):
    # counts_freq=1400
    # data_freq=325
    # Si=-0.7
    # S = S * (counts_freq / data_freq) ** Si
    a = np.array([a0, a1, a2, a3, a4, a5, a6, a7])
    ivals = np.arange(len(a))
    vals = np.zeros_like(S)
    logS = np.log10(S)
    for i in range(len(S)):
        vals[i] = np.dot(a, logS[i]**ivals)
    return vals

def tot_num_sources(Svals, survey_area):
    dNdSvals = 10 ** (SCs_Bondi(Svals)) * (Svals / 1000) ** (-2.5)
    integ = trapezoid(dNdSvals, Svals / 1000)
    number_sources = int(np.round(integ * (survey_area * (np.pi / 180) ** 2)))
    return number_sources

# Define the function to find the root for
def source_diff(max_flux, min_flux, target_num_sources, survey_area):
    Svals = np.linspace(min_flux, max_flux, 10000)
    num_sources = tot_num_sources(Svals, survey_area)
    return num_sources - target_num_sources

# Wrapper function to find max_flux
def find_max_flux(min_flux, target_num_sources, survey_area, initial_guess=10):
    # Check the initial values to ensure the root lies within the range
    lower_value = source_diff(min_flux, min_flux, target_num_sources, survey_area)
    upper_value = source_diff(initial_guess, min_flux, target_num_sources, survey_area)
    
    # Expand the search range if initial guess doesn't have different signs
    if lower_value * upper_value > 0:
        factor = 10  # Factor to expand the upper limit
        while lower_value * upper_value > 0:
            initial_guess *= factor
            upper_value = source_diff(initial_guess, min_flux, target_num_sources, survey_area)
            # Safety condition to avoid infinite loop
            if initial_guess > 1e2:  # Adjust this value based on expected max_flux range
                raise ValueError("Cannot find suitable bracket. Try adjusting min_flux or initial_guess.")

    # Use brentq to find the root (max_flux) within the bracket [min_flux, initial_guess]
    max_flux = brentq(
        source_diff,          # Function to find the root of
        min_flux,             # Lower bound of the search interval
        initial_guess,        # Upper bound of the search interval
        args=(min_flux, target_num_sources, survey_area)  # Additional arguments to pass to `source_diff`
    )
    return max_flux

# Example usage (replace with your specific function and values)
min_flux = 50e-3  # Example lower limit flux
survey_area = 1.26**2  # Example survey area in square degrees
target_num_sources = 2000  # Desired total number of sources

# Calculate the upper flux limit using `find_max_flux`
try:
    max_flux = find_max_flux(min_flux, target_num_sources, survey_area)
    print(f"Upper flux limit (max_flux) is: {max_flux:.6e}")
except ValueError as e:
    print(e)