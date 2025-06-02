# extra
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, trapezoid

def SCs_Mandal(S,  a0=1.655, a1=-0.1150, a2=0.2272, a3=0.51788, a4=-0.449661, a5=0.160265, a6=-0.028541, a7=0.002041):
    
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
    dNdSvals = (10 ** (SCs_Bondi(Svals))) * (Svals/1000)**(-2.5)
    integ = trapezoid(dNdSvals, Svals/1000)
    number_sources= int(np.round(integ*(survey_area*(np.pi/180)**2)))
    return(number_sources)

S_min=20e-3
S_max=10
S = np.linspace(1e-3, S_max, 10000)
Svals=np.linspace(S_min, S_max, 10000)
survey_area=1.5
number_sources=  tot_num_sources(Svals, survey_area)
print('number of sources', number_sources)


# plt.plot(S, 10**(SCs_Mandal(S)),label=r' Mandal 2021 $dN/dS \times S^{2.5}$')
#plt.plot(S, 10** (SCs_Mandal(S)) * (S/1000)**(-2.5),label='Mandal 2021 dN/dS')
#plt.plot(S, 10**(SCs_Bondi(S)),label=r' Bondi 2008 $dN/dS \times S^{2.5}$')
plt.plot(S, 10**SCs_Bondi(S) * (S/1000)**(-2.5),label=' Bondi 2008 dN/dS')

plt.xscale('log')
plt.yscale('log')
plt.axvline(x=20e-3, color='red', linestyle='--', label=rf'{S_min}')
#plt.axvline(x=50e-3, color='green', linestyle='--', label=r'50e-6')
plt.axvline(x=S_max, color='blue', linestyle='--', label=rf'{S_max}')
plt.xlabel('Flux S [mJy]',size=12)
plt.ylabel(' dN/dS',size=12)
plt.title('Source Counts')
plt.legend()
plt.show()
plt.close()
