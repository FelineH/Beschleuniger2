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



##############################################################################################################################################################################

#3.a)

#define functions

#universal function that should be used for the calculation later for the spectral density
def universal_function(omega):
    integral = integrate.quad(lambda x: kv(v_bessel,x), omega/omega_critical, np.inf) #https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html
    if integral[1] >= integral[0]*1e-3: # second entry holds numerical integration error
        raise Warning("Error of integration in range of value. Occured ad omega={:8f}".format(omega))
    S = (9*np.sqrt(3)/(8*np.pi))*(omega/omega_critical)*integral[0]
    return S



#calculation of spectral density 

y_values = [] # empty list --> plot
x_values = []
for i in np.arange(0.001,5.000,0.001): # range()-function can't be used to generate the range of float numbers --> arange() can have float numbers as inputs
    spectral_density = (P_s*universal_function(i*omega_critical))/omega_critical
    y_values.append(spectral_density) # add x and y-values to the empty list in order to plot the spectral density in the range from 0,001*omega_critical to 5*omega_critical
    x_values.append(i*omega_critical)    



#draw function

#draw angular spectral density in a linear form
fig, ax = plt.subplots()
ax.plot(x_values,y_values)
ax.axvspan(*minmax, color='C3', alpha=1  ,label='Empfindlichkeit der Diode')
ax.set_title('Spektrale Leistungsdichte')
ax.set_xlabel(r"$\omega$ in $10^{19} \, 1/s$")
ax.set_ylabel(r"$\mathrm{d}P/\mathrm{d}\omega$")
ax.grid(True, alpha = 0.7)
fig.savefig('build/spektrum_synchrotronstrahlung.pdf')
fig.clf()

#draw angular spectral density in a logarithmic form 
fig,ax = plt.subplots()
ax.plot(x_values,y_values)
ax.axvspan(*minmax, color='C3', alpha=0.3, label='Empfindlichkeit der Diode')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_title('Spektrale_Leistungsdichte doppelt-logarithmisch')
ax.set_xlabel(r"$\omega$ in $1/s$")
ax.set_ylabel(r"$\mathrm{d}P/\mathrm{d}\omega$")
ax.set_title('Spektrale Leistungsdichte')
ax.grid(True, alpha = 0.7)
fig.savefig('build/spektrum_synchrotronstrahlung_logarithmisch.pdf')
fig.clf()



#############################################################################################################################################################################

# 3.b)

spectral_density = np.array(y_values)  
frequency = np.array(x_values)
minmax = np.array([200e-9, 800e-9]) # spectral intervall of light which can pass through the vacuum window and the air, and can be measured at the photodiode
minmax = c/minmax * 2 * np.pi # wavelength to angular frequency


#Calculation of the total power

P_ges=np.sum(spectral_density) # ratio between power emitted and power measured is the same as the sum of all power per spectral intervall
P_diode=np.sum(spectral_density[np.logical_and(frequency < minmax[0], frequency > minmax[1])]) # proportion of the spectrum measured by the diode in ratio to the power per spectral 
                                                                                               # intervall between 200nm and 800nm


#Calculation of the power ratio P_ges at the photodiode

#due to geometry just a small ratio gets to the diode (analogous to  Blatt2 A2)
alpha = 0.01/10 # with sin alpha = alpha = a/L = 0.001; a = 1cm = 0.01m; L = 10m 
P_meas = P_diode/P_ges / (2 * np.pi) * alpha

print('Ratio:', P_meas)
