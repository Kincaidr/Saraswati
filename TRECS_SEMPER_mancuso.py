from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from SEMPER import SEMPER_SFG_AGN_counts
from TRECS_counts import source_counts
from Mancuso_2017_counts import Mancuso_counts


S_M,M_counts=Mancuso_counts()
Semper_M,Semper_counts=  SEMPER_SFG_AGN_counts()
TRECS_M,TRECS_count=  source_counts()
plt.plot(S_M,M_counts,label='Mancuso 2017 SFG + AGN model (Mancuso +2017)',color="#BB1616EB")
plt.plot(Semper_M,Semper_counts,label='SEMPER SFG model (SEMPER +2024)',color='green')
plt.scatter(TRECS_M,TRECS_count,label='TRECS SFG model (TRECS +2024)',color='orange')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()