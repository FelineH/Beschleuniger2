import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 
import scipy.integrate

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#for K in np.arange(1,2.5,0.5):
    #größen aus der aufgabenstellung
perioden = 20
λu = 0.1 #m
E_e = 1.5e9 #eV
ħω = np.linspace(1, 450, 450) #eV Photonenenergie E_Photon
K = 1.5

    #aus den parametern abgeleitete größen
ωu = (2 * np.pi / λu) * const.c # Formel aus Skript 
ω  = ħω * const.e / const.hbar # ω aus E_photon = hbar*ω bestimmen 
                                # ω ist Spaltenvektor mit 450 Einträgen (450,)
    #print('omega', ω, ω.shape)
γ = E_e/511e3

    #mittlere geschwindigkeit berechnen
β_m = 1- 1/(2 * γ**2) * (1 + K **2 / 2)

    #integrationsbereich also die "zeitachse" erstellen
Tu = perioden * λu / β_m / const.c
t = np.linspace(-Tu/2, Tu/2, 1000)

    #abstand zwischen ladung und beobachtungspunkt
r = - β_m * const.c * t + const.c * K **2 / (8*ωu*γ**2) * np.sin(2*ωu*t) #als Annahme aus Aufgabenstellung gegeben, wobei exp(i*r_p*ω) nur ein Phasenfaktor ist, 
                                                                            #der wegfällt mit r_p=0

    #vektoren aus der aufgabenstellung implementieren
n_vec = np.array((0,0,1)) #n_vec ist noch ein Zeilenvektor, der in für a, I_real und I_im transponiert werden muss
β_vec = np.array((K/γ * np.sin(ωu*t), 0*t, β_m + K**2/(4*γ**2) * np.cos(2*ωu*t))) #0*t damit die y-komp die gleiche dim hat wie x und z-komp
                                                                                    # beta_vec hat 3 Zeilen für x-,y- und z-Komponente mit je 1000 Einträgen, 
                                                                                    #d.h. Matrix der Form (3,1000) 
##################################################################################################################################################################
#d)
#β_vec = np.array((K/γ * np.sin(ωu*t), 0*t, np.full_like(t,β_m))) #longitudinale Komponente von β_vec ist in allen Zeitschritten konstant 
#print(β_vec, β_vec.shape)
#print(n_vec.shape)
##################################################################################################################################################################

    # integrand berechnen, in real und imaginärteil aufgeteilt
a = np.cross(n_vec.T, β_vec.T)
I_real = np.cross(n_vec.T, np.cross(n_vec.T, β_vec.T))[:,0] * np.cos(- np.outer(ω, (t + r / const.c))) #Funktion np.outer() berechnet äußeres Produkt der 
                                                                                                        #Vektoren ω und t + r / const.c
                                                                                                        # mit (t + r / const.c) als Spaltenvektoren mit 1000 Einträgen (1000,)
I_im = np.cross(n_vec.T, np.cross(n_vec.T, β_vec.T))[:,0] * np.sin(- np.outer(ω, (t + r / const.c)))

#E-Feld berechnen
#integrate.simps ist eine numerische methode der integration
E_real = scipy.integrate.simps(I_real, t)
E_im = scipy.integrate.simps(I_im, t)

#Intensität berechnen bzw. etwas das dazu proportional ist
A = (E_real+E_im)**2

#plotten
fig, ax = plt.subplots()
ax.plot(ω , A , )
ax.set_title('Spektrum')
ax.set_ylabel(r"$ E² \, \mathrm{in} \, V²s²/m²$")
ax.set_xlabel(r"$ ω \, \mathrm{in} \, s$")
    #ax.set_xlim(-1,1)
    #ax.set_ylim(-1,1)
#fig.savefig('spektrum-{:06.1f}.pdf'.format(K))
fig.savefig('build/spektrum_{}{}_{}{}.pdf'.format("K", K, "Perioden", perioden))


