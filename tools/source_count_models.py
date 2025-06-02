import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

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

def SCs_Bondi(S,  a0=0.805, a1=0.493, a2=0.564, a3=-0.129, a4=-0.195, a5=0.110, a6=-0.017):
    a=np.array([a0, a1, a2, a3, a4, a5, a6])
    ivals = np.arange(len(a))
    logS = np.log10(S)
    vals = np.zeros_like(S)
    for i in range(len(S)):
         vals[i] = np.dot(a, logS[i]**ivals)
    return vals


def tot_num_sources(Svals, survey_area):
    dNdSvals = 10 ** (SCs_Bondi(Svals)) * (Svals/1000)**(-2.5)
    integ = trapezoid(dNdSvals, Svals/1000)
    number_sources= int(np.round(integ*(survey_area*(np.pi/180)**2)))
    return(number_sources)

name='COSMOS'
min_flux=2.4*1e-3
max_flux=100
length=0.55
survey_area=length**2
Svals = np.linspace(min_flux, max_flux, 10000)
number_sources= tot_num_sources(Svals, survey_area)
print('Total number of sources for '+str(name)+ ' is '+ str(number_sources))