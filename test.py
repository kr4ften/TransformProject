import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

R = 5000
C = 20 * 10**(-9)
L = 50 * 10**(-3)

numerator = [C*L*R, 0, R]
denominator = [C*L*R, L, R]

LTI = signal.TransferFunction(numerator,denominator)
tDataLTI, yDataLTI = signal.step(LTI)

# Analytical expression
a = - 1 / (2*R*C)
b = np.sqrt(np.abs(L-4*R ** 2 * C /(4* R **2 * C ** 2 * L)))
c = np.sqrt(np.abs(L ** 2 - 4 * R ** 2 * C * L))

calcResponse = lambda t : 1 - 2*L / (c) * np.exp(t*a) * np.sin(b * t)

# Testing...
b2 = np.sqrt(1/(C*L) - 1/(4* R**2 * C**2))
c2 = np.sqrt(4 * R * C * L - L**2)

calcResponse2 = lambda t : 1 - 2*L/(c2) * np.exp(t*a) * np.sin(b2 * t)

# Plot function

tInput = np.linspace(min(tDataLTI),max(tDataLTI),1000)

plt.figure(figsize=(9, 3), dpi=100)
plt.plot(tDataLTI, yDataLTI, label='Numerical')
plt.plot(tInput, calcResponse(tInput), color='r', label='Analytical')
plt.legend()
plt.xlabel('Time [s]]')
plt.ylabel('Voltage [V]')
#plt.grid()
plt.show()