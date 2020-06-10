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
R_56 = 1000 #matrix element for 3c) normed by λ


#define functions

#electrondensity
def rho(elektronen,resolution = 500):
    x = np.linspace(elektronen[0].min(), elektronen[0].max(), resolution)
    dichte = np.zeros(resolution)
    for i in range(resolution-1):
        dichte[i] = np.sum((elektronen[0]<x[i+1])*(elektronen[0]>x[i]))
    return np.array((x, dichte/np.max(dichte))) # shape (2, 500)

#bunching factor
#def bunching(elektronen): 
#    bunching_h = []
#    bunching = 0
#    #harmonic = np.linspace(1,51,50) #array with 50 inputs 
#    #compute the bunchingfactor for the harmonics 1 to 50
#    for h in range (1,51,1):   #from 1 to 50
#        for i in range(1,3001,1): #from 1 to 3000
#            sum_function = np.exp(1j*2*np.pi*h*elektronen[0,i-1]) #start with x_data_uniform[i-1] because its an array with indexes form 0 to 2999 with electrons[0,i] = x_data_uniform[i-1]
#            bunching += sum_function 
#        bunching_factor = 1/N * np.absolute(bunching) #abs() for absolute value of the exponential function
#        bunching_h.append(bunching_factor)
#    return bunching_h

def bunching(elektronen): 
    bunching_h = []
    bunching = 0
    for h in range (1,51,1):
        for i in range(1,3001,1):
            bunching = 1/N*abs(np.sum(np.exp(1.j*2*np.pi*h*elektronen[0,i-1])))
        bunching_h.append(bunching)
    return bunching_h

##bunching factor
#bunching_h = []
#bunching = 0
#harmonic = np.linspace(1,51,50) #array with 50 inputs 

#compute the bunchingfactor for the harmonics 1 to 50
#for h in range (1,51,1):   #from 1 to 50
 #   for i in range(1,3001,1): #from 1 to 3000
  #      sum_function = np.exp(1j*2*np.pi*h*x_data_uniform[i-1]) #start with x_data_uniform[i-1] because its an array with indexes form 0 to 2999
   #     #print(sum_function)
    #    bunching += sum_function 
     #   #print(bunching)
   # bunching_factor = 1/N * np.absolute(bunching) #abs() for absolute value of the exponential function
   # bunching_h.append(bunching_factor)


#a)
########################################################################################################################################################################
#energy deviation
x_data_uniform = np.random.uniform(0,2,N) #random homogenous distribution of electron in longitudinal direction
y_data_normal = norm.rvs(size =N, loc=0, scale=0.001) 
        #y_data_normal= np.random.normal(0, 0.001, N)
        #generate normally distributed random variable using scipy.stats module's norms.rvs(); loc argument corresponds to the meean
        #and scale to standard deviation and size to the number of random variates

electrons = np.array([x_data_uniform, y_data_normal]) # shape (2,3000)
#electrons[0,:]) is the same as (x_data_uniform)

#electrondensitiy for a)
dichte = rho(electrons) # to set in plot function

#bunchingfactor for a)
bunchingfactor = bunching(electrons)
harmonic = np.linspace(1,51,50) #array with 50 inputs
print(bunchingfactor)

#plot DeltaE/E
fig, ax = plt.subplots(nrows=3, ncols=1,figsize=(15,15))
ax[0].scatter(x_data_uniform,y_data_normal, marker='.', color="red")
ax[0].set_ylabel(r"$ \Delta E/E$")
ax[0].set_xlabel(r"$ z/\lambda$")
ax[0].set_xlim(0,2)
ax[0].set_ylim(-0.01,0.01)
#fig.savefig('build/normalverteilte_Elektronen_a.pdf')

#plot Elektronendichtelectron density
#fig, ax = plt.subplots(1, 1)
#ax.plot(x_data_uniform, density_electron, '-', color="red")
ax[1].plot(dichte[0], dichte[1], color="red")
ax[1].set_ylabel(r"$ \rho(z/\lambda)$")
ax[1].set_xlabel(r"$ z/\lambda$")
ax[1].set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
#fig.savefig('build/electron_density_a.pdf')

#plot bunching factor 
#fig, ax = plt.subplots(1, 1)
ax[2].plot(harmonic, bunchingfactor, 'o-', color='red')
ax[2].set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax[2].set_xlabel(r"$\mathrm{Harmonische} \, h$")
#fig.savefig('build/bunching_factor_a.pdf')

plt.title('Zufällige Energieverteilung')
plt.tight_layout()
plt.savefig('build/a.pdf')


#b)
########################################################################################################################################################################
#energy deviation
y_data_normal_sinus = y_data_normal + A*np.sin(2*np.pi*x_data_uniform) #energy of electrons is sinusoidally distributed

electrons_sinus = np.array([x_data_uniform, y_data_normal_sinus]) # shape (2,3000)

#######################################################################
#no changes from task a to b for electron density and bunching factor
#######################################################################
#electrondensity for b)
dichte_b = rho(electrons_sinus) 

#bunchingfactor for b)
bunchingfactor_b = bunching(electrons_sinus)
harmonic = np.linspace(1,51,50) #array with 50 inputs

#plot DeltaE/E
fig, ax = plt.subplots(nrows=3, ncols=1,figsize=(15,15))
ax[0].scatter(x_data_uniform,y_data_normal_sinus, marker='.', color="red")
ax[0].set_ylabel(r"$ \Delta E/E$")
ax[0].set_xlabel(r"$ z/\lambda $")
ax[0].set_xlim(0,2)
ax[0].set_ylim(-0.01,0.01)
#fig.savefig('build/normalverteilte_Elektronen_b.pdf')

#plot electron density
#fig, ax = plt.subplots(1, 1)
ax[1].plot(dichte_b[0], dichte_b[1], '-', color="red")
ax[1].set_ylabel(r"$ \rho(z/\lambda)$")
ax[1].set_xlabel(r"$ z/\lambda$")
ax[1].set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
#fig.savefig('build/electron_density_b.pdf')

#plot bunching factor 
#fig, ax = plt.subplots(1, 1)
ax[2].plot(harmonic, bunchingfactor_b, 'o-', color='red')
ax[2].set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax[2].set_xlabel(r"$\mathrm{Harmonische} \, h$")
#fig.savefig('build/bunching_factor_b.pdf')

plt.title('Energiemodulation')
plt.tight_layout()
plt.savefig('build/b.pdf')

#c)/d)
########################################################################################################################################################################
#chicane
x_data_uniform_chicane = (x_data_uniform + y_data_normal_sinus * R_56) #path length of electrons is dependent from the energy deviation from b) (y_data_normal_sinus) 
    	                                                                #and from the transfer matrix element R_56

electrons_chicane = np.array([x_data_uniform_chicane, y_data_normal_sinus]) # shape (2,3000)

#electrondensity for c/d)
dichte_cd = rho(electrons_chicane) 

#bunchingfactor for c/d)
bunchingfactor_cd = bunching(electrons_chicane)
harmonic = np.linspace(1,51,50) #array with 50 inputs

#########################################################################
#########################################################################

#plot DeltaE/E
fig, ax = plt.subplots(nrows=3, ncols=1,figsize=(15,15))
ax[0].scatter(x_data_uniform_chicane,y_data_normal_sinus, marker='.', color="red")
ax[0].set_ylabel(r"$ \Delta E/E $")
ax[0].set_xlabel(r"$ z/\lambda$")
ax[0].set_xlim(0,2)
ax[0].set_ylim(-0.01,0.01)
#fig.savefig('build/normalverteilte_Elektronen_d.pdf')

#plot electron density
#fig, ax = plt.subplots(1, 1)
ax[1].plot(dichte_cd[0], dichte_cd[1], '-', color="red")
ax[1].set_ylabel(r"$ \rho(z/\lambda)$")
ax[1].set_xlabel(r"$ z/\lambda $")
ax[1].set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
#fig.savefig('build/electron_density_d.pdf')

#plot bunching factor 
#fig, ax = plt.subplots(1, 1)
ax[2].plot(harmonic, bunchingfactor_cd, 'o-', color='red')
ax[2].set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax[2].set_xlabel(r"$\mathrm{Harmonische} \, h$")
#fig.savefig('build/bunching_factor_d.pdf')

plt.title('Schikane und Dichtemodulation')
plt.tight_layout()
plt.savefig('build/cd.pdf')

#e)/f)
########################################################################################################################################################################
#energy deviation
y_data_normal_sinus_sinus = y_data_normal_sinus + A*np.sin(2*np.pi*x_data_uniform_chicane) # modulation of electrons' energy is sinusoidally distributed

x_data_uniform_chicane_new = (x_data_uniform + y_data_normal_sinus * 25) # with new R_56 value (25) to cause a density modulation with high harmonics of the laser's wave length
electrons_sinus_sinus = np.array([x_data_uniform_chicane_new, y_data_normal_sinus_sinus])

#electrondensity for e/f)
dichte_ef = rho(electrons_sinus_sinus) 

#bunching factor for e/f)
bunchingfactor_ef = bunching(electrons_sinus_sinus)
harmonic = np.linspace(1,51,50) #array with 50 inputs

#########################################################################
#########################################################################

#plot DeltaE/E
fig, ax = plt.subplots(nrows=3, ncols=1,figsize=(15,15))
ax[0].scatter(x_data_uniform_chicane,y_data_normal_sinus_sinus, marker='.', color="red")
ax[0].set_ylabel(r"$ \Delta E/E$")
ax[0].set_xlabel(r"$ z/\lambda $")
ax[0].set_xlim(0,2)
ax[0].set_ylim(-0.01,0.01)
#fig.savefig('build/normalverteilte_Elektronen_fR1000.pdf')

#plot Elektronendichtelectron density
#fig, ax = plt.subplots(1, 1)
ax[1].plot(dichte_ef[0], dichte_ef[1], '-', color="red")
ax[1].set_ylabel(r"$ \rho(z/\lambda) $")
ax[1].set_xlabel(r"$ z/\lambda $")
ax[1].set_xlim(0,2)
#ax.set_ylim(-0.01,0.01)
#fig.savefig('build/electron_density_fR1000.pdf')

#plot bunching factor 
#fig, ax = plt.subplots(1, 1)
ax[2].plot(harmonic, bunchingfactor_ef, 'o-', color='red')
ax[2].set_ylabel(r"$\mathrm{Bunching-Faktor} \, b_{h}$")
ax[2].set_xlabel(r"$\mathrm{Harmonische} \, h$")
#fig.savefig('build/bunching_factor_fR1000.pdf')

plt.title('Zwei Schikanen und Dichtemodulation')
plt.tight_layout()
plt.savefig('build/ef.pdf')