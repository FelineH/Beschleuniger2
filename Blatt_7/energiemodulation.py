import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 
from scipy.stats import norm

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#define constants
N = 3000 #number of electrons
A = 0.003 #amplitude for 3b) 

#a)
########################################################################################################################################################################
#energy deviation
x_data_uniform = np.random.uniform(0,2,N) #random homogenous distribution of electron in longitudinal direction
y_data_normal = norm.rvs(size =N, loc=0, scale=0.001) 
#generate normally distributed random variable using scipy.stats module's norms.rvs(); loc argument corresponds to the meean
#and scale to standard deviation and size to the number of random variates


#electrondensity
density_electron = x_data_uniform/N 


#bunching factor
bunching_h = []
bunching = 0
harmonic = np.linspace(1,51,50) #array with 50 inputs 

#compute the bunchingfactor for the harmonics 1 to 50
for h in range (1,51,1):   #from 1 to 50
    for i in range(1,3001,1): #from 1 to 3000
        sum_function = np.exp(1j*2*np.pi*h*x_data_uniform[i-1]) #start with x_data_uniform[i-1] because its an array with indexes form 0 to 2999
        #print(sum_function)
        bunching += sum_function 
        #print(bunching)
    bunching_factor = 1/N * np.absolute(bunching) #abs() for absolute value of the exponential function
    bunching_h.append(bunching_factor)
#print(bunching_h)

#plot DeltaE/E
fig, ax = plt.subplots(1, 1)
plt.scatter(x_data_uniform,y_data_normal, marker='.', color="red")
ax.set_ylabel(r"$ \Delta E/E \, \mathrm{in} \, \%$")
ax.set_xlabel(r"$ z/\lambda \, \mathrm{in} \, m$")
ax.set_xlim(0,2)
ax.set_ylim(-0.01,0.01)
fig.savefig('build/normalverteilte_Elektronen_a.pdf')

#plot Elektronendichtelectron density
fig, ax = plt.subplots(1, 1)
ax.plot(x_data_uniform, density_electron, '-', color="red")
ax.set_ylabel(r"$ \rho(z/\lambda)$")
ax.set_xlabel(r"$ z/\lambda \, \mathrm{in} \, m$")
ax.set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
fig.savefig('build/electron_density_a.pdf')

#plot bunching factor 
fig, ax = plt.subplots(1, 1)
ax.plot(harmonic, bunching_h, 'o-', color='red')
ax.set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax.set_xlabel(r"$\mathrm{Harmonische} \, h$")
fig.savefig('build/bunching_factor_a.pdf')

#b)
########################################################################################################################################################################
#energy deviation
y_data_normal_sinus = y_data_normal + A*np.sin(2*np.pi*x_data_uniform) #energy of electrons is sinusoidally distributed


#######################################################################
#no changes from task a to b for electron density and bunching factor
#######################################################################
#electrondensity
density_electron_sinus = x_data_uniform/N 

#bunching factor
bunching_h = []
bunching = 0
harmonic = np.linspace(1,51,50) #array with 50 inputs 

#compute the bunchingfactor for the harmonics 1 to 50
for h in range (1,51,1):   #from 1 to 50
    for i in range(1,3001,1): #from 1 to 3000
        sum_function = np.exp(1j*2*np.pi*h*x_data_uniform[i-1]) #start with x_data_uniform[i-1] because its an array with indexes form 0 to 2999
        #print(sum_function)
        bunching += sum_function 
        #print(bunching)
    bunching_factor = 1/N * np.absolute(bunching) #abs() for absolute value of the exponential function
    bunching_h.append(bunching_factor)
#print(bunching_h)

#########################################################################
#########################################################################

#plot DeltaE/E
fig, ax = plt.subplots(1, 1)
plt.scatter(x_data_uniform,y_data_normal_sinus, marker='.', color="red")
ax.set_ylabel(r"$ \Delta E/E \, \mathrm{in} \, \%$")
ax.set_xlabel(r"$ z/\lambda \, \mathrm{in} \, m$")
ax.set_xlim(0,2)
ax.set_ylim(-0.01,0.01)
fig.savefig('build/normalverteilte_Elektronen_b.pdf')

#plot Elektronendichtelectron density
fig, ax = plt.subplots(1, 1)
ax.plot(x_data_uniform, density_electron_sinus, '-', color="red")
ax.set_ylabel(r"$ \rho(z/\lambda)$")
ax.set_xlabel(r"$ z/\lambda \, \mathrm{in} \, m$")
ax.set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
fig.savefig('build/electron_density_b.pdf')

#plot bunching factor 
fig, ax = plt.subplots(1, 1)
ax.plot(harmonic, bunching_h, 'o-', color='red')
ax.set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax.set_xlabel(r"$\mathrm{Harmonische} \, h$")
fig.savefig('build/bunching_factor_b.pdf')