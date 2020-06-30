import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as const

T0 = 384e-9 #s
I = 10e-3 #A
E = 1.5e9 #eV
sigma_t = 40e-12 #s
r = 0.01 #m
sigma = 1.4e6 #Ohm⁻¹m⁻¹

c = const.c
e = const.e
epsilon_0 = const.epsilon_0

omega_0 = 2 * np.pi / T0
N = I * T0 / e
L = c * T0

rate = np.zeros(180)
    
for index, arbeitspunkt in enumerate(np.linspace(3.05, 3.95, 180)):
    omega_b = omega_0 * arbeitspunkt

    p = np.round(np.linspace(-1000, 1000, 2000))

    ReZ = np.sign(omega_0*p+omega_b) * L/(np.pi*r**3 * np.sqrt(2* epsilon_0 * sigma * np.absolute(omega_0*p+omega_b)))

    F = np.exp(-0.5*(omega_0*p+omega_b)**2*sigma_t**2)

    summe = np.sum(ReZ*F)

    rate[index] = -N*e**2*c/(2*omega_b*E*T0**2)*summe


#fig3,ax = plt.subplots()
#ax.plot(freq, signal_fft.real)
#ax.set_title(r"FFT-Signal nach 1000 Perioden")
#ax.set_xlabel(r"ω in $1/ns$")
#ax.set_ylabel(r"FFT-Signal")
#fig3.savefig(r"Signal_Zeit_FFT.pdf")
plt.plot(np.linspace(3.05, 3.95, 180), rate)
plt.ylabel("Anstiegs-/Abfallrate s")
plt.xlabel("Arbeitspunkt")
plt.savefig("Arbeitspunkt.pdf")