import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from astropy.cosmology import Planck18 as cosmo
from scipy.interpolate import interp1d
# Given parameters for the luminosity function
import numpy as np
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# Constants
H0 = 70
Om0 = 0.3
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

def FSRQ_cont():
    a = 0.760
    b = 2.508
    n_0 = (10**(-10.382))*u.Mpc**-3 *u.dex**-1  # Mpc^-3
    n_0=n_0.to(u.m**-3*u.dex**-1 )
    l_0 = 10**34.323* u.W / u.Hz
    k_evo = -0.996
    z_top_0 = 1.882
    delta_z_top = 0.018
    m_ev = -0.166
    return a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top

def BLLac_cont():
    a = 0.723
    b = 1.618
    n_0 = (10**(-6.879))*u.Mpc**-3*u.dex**-1 # Mpc^-3
    n_0=n_0.to(u.m**-3*u.dex**-1 )
    l_0 = 10**32.638* u.W / u.Hz
    k_evo = 0.208
    z_top_0 = 1.282
    delta_z_top = 1
    m_ev = 1
    return a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top

def SS_AGNs_cont():
    a = 0.559
    b = 2.261
    n_0 = (10**(-5.970))*u.Mpc**-3*u.dex**-1   # Mpc^-3
    n_0=n_0.to(u.m**-3*u.dex**-1 )
    l_0 = 10**32.490* u.W / u.Hz
    k_evo = 1.226
    z_top_0 = 0.977
    delta_z_top = 0.842
    m_ev = 0.282
    return a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top

def z_top(L,z_top_0, delta_z_top, l_0):
    L = L.to(l_0.unit)
    return z_top_0 + delta_z_top / (1 + l_0 / L)

def L_star(z, L,k_evo, l_0, m_ev):
    zt = z_top(L,z_top_0,delta_z_top,l_0)
    exponent = k_evo * z * (2 * zt - 2 * (z**m_ev * zt**(1 - m_ev))) / (1 + m_ev)
    return l_0 * 10**(exponent)

def phi_L(L, z,a, b, n_0):
    Ls = L_star(z, L,k_evo,l_0,m_ev)
    return n_0 / ((L / Ls)**a + (L / Ls)**b)

def integrand(z, S):
    S=S*u.Jy
    S_si = S.to(u.W / (u.m**2 * u.Hz))
    d_L = cosmo.luminosity_distance(z).to(u.m) 
    L = 4 * np.pi * d_L**2 * S_si #* (1 + z)**(1 + alpha) # in W/Hz
    #logL= np.log10(L)
    dV_dz = cosmo.differential_comoving_volume(z).to(u.m**3 / u.sr) # m^3 / sr
    area=(6.5*u.deg**2).to(u.sr)
    return (phi_L(L, z,a,b,n_0) * (1/(S*np.log(10))) *area* dV_dz).value

def dNdS(S):
    result, _ = quad(integrand, 0.001, 2, args=(S,), limit=200, epsabs=1e-6, epsrel=1e-6)
    return (result)

def plot_luminosity_function(S_vals,counts):
    euclidean_norm = S_vals**2.5 * counts
    plt.figure()
    plt.plot(np.log10(S_vals), np.log10(euclidean_norm))
    #plt.xscale('log')
    #plt.ylim(1e-5, 1e3)
    plt.xlabel('Flux Density S [Jy]')
    plt.ylabel(r'$S^{2.5} \frac{dN}{dS}$ [sr$^{-1}$ Jy$^{1.5}$]')
    plt.grid(True, which="both", ls="--")
    plt.title('Euclidean-normalized Source Counts')
    plt.show()

if __name__ == "__main__":
    S_vals = np.logspace(-5, 3, 100) 
    a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top= FSRQ_cont()
    counts_FSRQ = np.array([dNdS(S) for S in S_vals])
    a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top= BLLac_cont()
    counts_BLLac = np.array([dNdS(S) for S in S_vals])
    a, b, n_0, l_0, k_evo, m_ev, z_top_0, delta_z_top= SS_AGNs_cont()
    counts_SS_AGNs = np.array([dNdS(S) for S in S_vals])
    counts= counts_FSRQ + counts_BLLac + counts_SS_AGNs
    plot_luminosity_function(S_vals,counts)