import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 
import scipy.integrate

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#define Constants

λu = 0.1 #meter
E_e = 1.5e9 #eV
ħω = np.linspace(1, 450, 450) #eV #photon energy E_Photon
#ϑ  = np.linspace(0, 2*np.pi/10000, 50)
#φ  = np.linspace(0, 2*np.pi/10000, 50)
#ϑ = 0.9
#φ = np.pi/2
d = 1 #Biegeradius
#n_vec = np.array((0,0,1)) #n_vec is a row vector, which has to be transposed for I_real und I_im 



#computed constants

ωu = (2 * np.pi / λu) * const.c # second  #Formel aus Skript 
ω  = ħω * const.e / const.hbar # second #ω computed by E_photon = hbar*ω
                                # ω is a column vector with 450 entries 
                                #shape (450,)
γ = E_e/511e3




#define functions

def beta_averaged(K):
    '''mittlere Geschwindigkeit berechnen in meter/second'''
    β_m = 1- 1/(2 * γ**2) * (1 + K **2 / 2)
    return β_m

def integration_interval(periods,K):
    '''Integrationsbereich bzw. die "Zeitachse" erstellen'''
    Tu = periods * λu / beta_averaged(K) / const.c
    t = np.linspace(-Tu/2, Tu/2, 1000)
    return t

#def r_scal_old(K, periods):    
#    '''Abstand zwischen Ladung und Beobachtungspunkt in meter'''
#    r = - beta_averaged(K) * const.c * integration_interval(periods,K) + const.c * K **2 / (8*ωu*γ**2) * np.sin(2*ωu*integration_interval(periods,K)) 
#    #as the assumption of the task, whereby exp(i*r_p*ω) is a phase factor (drops out because of r_p=0 )
#    return r
##print("r_scal_old(1.5,20)", r_scal_old(1.5,20), r_scal_old(1.5,20).shape)                                                                  

def beta_vec(K, periods):
    '''Vektoren aus der Aufgabenstellung implementieren'''
    β_vec = np.array((K/γ * np.sin(ωu*integration_interval(periods, K)), 0*integration_interval(periods, K), beta_averaged(K) + K**2/(4*γ**2) * np.cos(2*ωu*integration_interval(periods, K)))) 
    #0*t to align the dimension of y-element with the dimension of the x- and z-elements
    #beta_vec hat 3 rows for the  x-,y- und z-elements, each with 1000 entries --> shape (3,1000) 
    return β_vec

def r_scal(K, periods, theta, phi):
    '''Abstand zwischen Ladung und Beobachtungspunkt in meter'''
    re = np.array((const.c * K /γ/ωu * np.cos(ωu*integration_interval(periods,K)), 
                    0*integration_interval(periods,K),
                    - beta_averaged(K) * const.c * integration_interval(periods,K) - const.c * K **2 / (8*ωu*γ**2) * np.sin(2*ωu*integration_interval(periods,K)) )) 
    #position vektor of the charge
    #re from integration of β_vec
    rp = np.array((d * np.sin(theta) * np.cos(phi), d * np.sin(theta) * np.sin(phi), d * np.cos(theta)))[:,np.newaxis] #position vector of the observation point in spherical coordinates
    #np.array()[:,newaxis] creates a column vector (3,1) out of 1D array of shape (3,)
    r_vek = re + rp
    r = np.linalg.norm(r_vek, axis=0) #axis = 0, because input is a tuple
    return r # r with shape(1000,)

def n_vec(K, periods, theta, phi):
    '''Normalenvektor'''
    re = np.array((-const.c * K /γ/ωu * np.cos(2*ωu*integration_interval(periods,K)), 0*integration_interval(periods,K),- beta_averaged(K) * const.c * integration_interval(periods,K) + const.c * K **2 / (8*ωu*γ**2) * np.sin(2*ωu*integration_interval(periods,K)) ))
    rp = np.array((d * np.sin(theta) * np.cos(phi), d * np.sin(theta) * np.sin(phi), d * np.cos(theta)))[:,np.newaxis]
    r_vek = re + rp
    n = r_vek/np.linalg.norm(r_vek, axis=0) # n = r/|r|
    return n


##################################################################################################################################################################
#d)
#β_vec = np.array((K/γ * np.sin(ωu*integration_interval(periods,K)), 0*integration_interval(periods,K), np.full_like(t,beta_averaged(K)))) 
#with a constant longitudinal component of β_vec in all time steps  
#print(β_vec, β_vec.shape)
#print(n_vec.shape)
##################################################################################################################################################################

#for loop to compute integrand, seperated in imaginary part and real part

φ = 0
for ϑ in np.linspace(0, 2*np.pi/1000, 50):
        #for K in range(1,2):
        #    for P in range(20,21):
        a = np.cross(n_vec(1, 20, ϑ, φ).T, beta_vec(1, 20).T)
        I_real = np.cross(n_vec(1, 20, ϑ, φ).T, a)[:,0] * np.cos(- np.outer(ω, (integration_interval(1, 20) + r_scal(1, 20, ϑ, φ) / const.c))) 
        #function np.outer() computes the outer product of vectors ω and (t + r / const.c) respectively of vectors ω and (integration_interval(P,K) + r_scal(K, P) / const.c)
        # (t + r / const.c) as a column vector with 1000 entries --> shape (1000,)
        I_im = np.cross(n_vec(1, 20, ϑ, φ).T, a)[:,0] * np.sin(- np.outer(ω, (integration_interval(1, 20) + r_scal(1, 20, ϑ, φ) / const.c)))

        #compute electrical field
        E_real = scipy.integrate.simps(I_real, integration_interval(1, 20)) # #integrate.simps is a numerical integration method
        E_im = scipy.integrate.simps(I_im, integration_interval(1, 20))

        #compute intensity 
        A = (E_real+E_im)**2 # E^2

        #compute omega
        K=1
        λ = λu/(2*γ**2)*(1+(K**2/2)+γ**2*ϑ**2) #fundamentale Undulatorlinie
        deltaλ = (λu*ϑ**2)/2 # winkelabhängige Rotverschiebung 
        ω_fundamental = (const.c*2*np.pi/λ)  
        abweichung_plus = (const.c*2*np.pi/(λ + deltaλ))  
        abweichung_minus = (const.c*2*np.pi/(λ - deltaλ))    


        #plot electrical field
        fig, ax = plt.subplots()
        ax.plot(ω , A , )
        plt.axvline(x = ω_fundamental, color='green', label ='Frequenz der fundamentalen Undulatorlinie')
        plt.axvline(x = abweichung_plus, color='r', label = 'Abweichung')
        plt.axvline(x = abweichung_minus, color='r', label = 'Abweichung')
        ax.set_title('Spektrum')
        ax.set_ylabel(r"$ E² \, \mathrm{in} \, V²s²/m²$")
        ax.set_xlabel(r"$ ω \, \mathrm{in} \, 1/s$")
        plt.legend()
        #ax.set_xlim(-1,1)
        #ax.set_ylim(-1,1)
        fig.savefig('build/spektrum_{}{}_{}{}_neu.pdf'.format("theta", ϑ, "phi", φ))