
from SEMPER import SEMPER_SFG_AGN_counts
import matplotlib.pyplot as plt
import numpy as np

Semper_M,Semper_counts=  SEMPER_SFG_AGN_counts()

x=np.log10(Semper_M)
y=np.log10(Semper_counts)

coeffs = np.polyfit(x, y, 8)
p = np.poly1d(coeffs)

print("Coefficients of the 7th order polynomial fit:", coeffs)
x_fit = np.linspace(min(x), max(x), 1000)
y_fit = p(x_fit)
plt.scatter(x, y, label='Data')
plt.plot(x_fit, y_fit, 'r-', label='7th Order Fit')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('7th Order Polynomial Fit')
plt.grid(True)
plt.show()