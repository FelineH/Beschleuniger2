import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 
import scipy.integrate

#größen aus der aufgabenstellung
perioden = 20
λu = 0.1 #m
E_e = 1.5e9 #eV
ħω = np.linspace(1, 450, 450) #eV
K = 1.5

#aus den parametern abgeleitete größen
ωu = 2 * np.pi / λu * const.c
ω  = ħω * const.e / const.hbar
γ = E_e/511e3

#mittlere geschwindigkeit berechnen
β_m = 1- 1/(2 * γ**2) * (1 + K **2 / 2)

#integrationsbereich also die "zeitachse" erstellen
Tu = perioden * λu / β_m / const.c
t = np.linspace(-Tu/2, Tu/2, 1000)

#abstand zwischen ladung und beobachtungspunkt
r = - β_m * const.c * t + const.c * K **2 / (8*ωu*γ**2) * np.sin(2*ωu*t)

#vektoren aus der aufgabenstellung implementieren
n_vek = np.array((0,0,1))
β_vek = np.array((K/γ * np.sin(ωu*t), 0*t, β_m + K**2/(4*γ**2) * np.cos(2*ωu*t))) #0*t damit die y-komp die gleiche dim hat wie x und z-komp

#integrand berechnen, in real und imaginärteil aufgeteilt
a = np.cross(n_vek.T, β_vek.T)
I_real = np.cross(n_vek.T, np.cross(n_vek.T, β_vek.T))[:,0] * np.cos(- np.outer(ω, (t + r / const.c)))
I_im = np.cross(n_vek.T, np.cross(n_vek.T, β_vek.T))[:,0] * np.sin(- np.outer(ω, (t + r / const.c)))

#E-Feld berechnen
#integrate.simps ist eine numerische methode der integration
E_real = scipy.integrate.simps(I_real, t)
E_im = scipy.integrate.simps(I_im, t)

#Intensität berechnen bzw. etwas das dazu proportional ist
A = (E_real+E_im)**2

#plotten
fig, ax = plt.subplots()
plt.plot(ω , A , )
ax.set_title('Spektrum')
ax.set_ylabel(r"$ E² /V²/m²$")
ax.set_xlabel(r"$ ω /s$")
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)
fig.savefig('spektrum.pdf')
#fig.clf()

