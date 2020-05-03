import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
from scipy.special import kv
import scipy.integrate as integrate

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#define constants
c = 299792458 #meter per second (speed of light)
e = 1.602176634*10**-19 #coulomb C (elementary charge)
epsilon_zero = 8.8541878128*10**-12 #As/Vm (vacuum permittivity)
m_electron=9.1093837015*10**-31#kg (mass of the electron)

R_bend = 3.33 #meter (Bending radius of thwe dipole magnet) 
E_electronbeam = 1.5 #GeV (energy of the electron beam)
I_beam = 100*10**-3 #Ampere A (beam current óf the emitted synchrotron light)
U_ring = 115.2 #meter (perimeter of the storage ring) 

v_bessel = 5/3 #order of the Bessel function and first paramter for th emodified Bessel function kv(v_bessel, x)




#computed constants
E_synchrotron = (88.5*E_electronbeam**4)/R_bend #keV with input E_electron in GeV (energy loss pro perimeter and electron)

gamma = (E_electronbeam*10**9)/(511*10**3)  #Lorentz factor with the electron's rest mass of 511 keV  #alternative calculation 1/(np.sqrt(1-(v/c)**2))

omega_critical = (3*c*(gamma**3))/(2*R_bend) # 1/second (critical angular frquency)
T_perimeter = U_ring/c #second (time for one perimeter in the storage ring)
#omega = (2*np.pi)/(T_perimeter) #1/second (angular frequency)

P_s = ((e**2)*c/(6*np.pi*epsilon_zero*(m_electron*c**2)**4))*((E_electronbeam*e)**4/(R_bend**2)) # watt W (angular integrated total power)

print("E_synchrotron =", E_synchrotron, "gamma =", gamma, "omega_crit = ", omega_critical,"angular integrated total power =", P_s)




#define functions

#universal function that should be used for the calculation later for the spectral density
def universal_function(omega):
    integral = integrate.quad(lambda x: kv(v_bessel,x), omega/omega_critical, np.inf)
    if integral[1] >= integral[0]*1e-3: # second entry holds numerical integration error
        raise Warning("Error of integration in range of value. Occured ad omega={:8f}".format(omega))
    S = (9*np.sqrt(3)/(8*np.pi))*(omega/omega_critical)*integral[0]
    return S


y_values = [] # empty list --> plot
x_values = []
for i in np.arange(0.001,5.000,0.001): # range()-function can't be used to generate the range of float numbers
    spectral_density = (P_s*universal_function(i*omega_critical))/omega_critical
    y_values.append(spectral_density)
    x_values.append(i*omega_critical)    

#proportion of the spectrum measured by the diode
#ratio between power emitted and power measured is the same as the sum of all power per spectral intervall
#in ratio to the power per spectral intervall between 200nm and 800nm

spectral_density = np.array(y_values)
frequency = np.array(x_values)
minmax = np.array([200e-9, 800e-9])
minmax = c/minmax * 2 * np.pi#wavelength to angular frequency
print(minmax)
#print(x_values)

P_ges=np.sum(spectral_density)
P_diode=np.sum(spectral_density[np.logical_and(frequency > minmax[0], frequency < minmax[1])])
print(P_ges)
print(P_diode)
print(P_diode/P_ges)


#draw function

#draw angular spectral density in a linear form
fig, ax = plt.subplots()
ax.plot(x_values,y_values)
ax.set_title('Spektrale Leistungsdichte')
ax.set_xlabel(r"$\omega$ in $10^{19} \, 1/s$")
ax.set_ylabel(r"$\mathrm{d}P/\mathrm{d}\omega$")
ax.grid(True, alpha = 0.7)
fig.savefig('build/spektrum_synchrotronstrahlung.pdf')
fig.clf()

#draw angular spectral density in a logarithmic form 
fig,ax = plt.subplots()
ax.plot(x_values,y_values)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title('Spektrale_Leistungsdichte doppelt-logarithmisch')
ax.set_xlabel(r"$\omega$ in $1/s$")
ax.set_ylabel(r"$\mathrm{d}P/\mathrm{d}\omega$")
ax.set_title('Spektrale Leistungsdichte')
ax.grid(True, alpha = 0.7)
fig.savefig('build/spektrum_synchrotronstrahlung_logarithmisch.pdf')
fig.clf()