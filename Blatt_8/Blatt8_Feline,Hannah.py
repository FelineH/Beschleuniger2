import numpy as np
import matplotlib.pyplot as plt

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#Konstanten
λ_u = 0.2 #m
E_0 = 1e9 #V/m
K = 2
ɣ = 3000
m_e = 9.1e-31 #kg
c = 3e8 #m/s
e = 1.6e-19 #C

N = 100 #anzahl elektronen
x = 0.01 #m, wegschritte

a1 = 4 * np.pi / λ_u
a2 = e * E_0 * K / (2 * m_e * c**2 * ɣ**2)

gain = np.zeros(21) #für c


#Differentialgleichungen
def dΨ(η, a1,):
    return a1 * η

def dη(Ψ, a2):
    return a2 * np.sin(Ψ)


#Startwerte definieren
Ψ = np.zeros((N, 1000))
Ψ[:,0] = np.random.uniform(0, 2*np.pi, N) #gleichverteilte Phase, von 0 bis 2pi, weil es schöner dargestellt wird, ändert nichts an berechnung
η = np.zeros((N, 1000))

for index, k in enumerate(np.linspace(-0.01, 0.01, 21)):
    k = np.round(k, 3) #durch nummerische ungenauigkeiten sind werte nicht genau die richtigen, deshalb runden
    η[:,0] = k  
#Runge-Kutta-Verfahren
    for i in range(1, 1000):

        j1 = x * dΨ(η[:,i-1], a1,)
        k1 = x * dη(Ψ[:,i-1], a2)
        j2 = x * dΨ(η[:,i-1]+k1/2, a1,)
        k2 = x * dη(Ψ[:,i-1]+j1/2, a2)
        j3 = x * dΨ(η[:,i-1]+k2/2, a1,)
        k3 = x * dη(Ψ[:,i-1]+j2/2, a2)
        j4 = x * dΨ(η[:,i-1]+k3, a1,)
        k4 = x * dη(Ψ[:,i-1]+j3, a2)

        Ψ[:,i] = Ψ[:,i-1] + 1/6 * (j1+ 2*j2 + 2*j3 + j4)
        Ψ[:,i] = Ψ[:,i]%(2*np.pi) #Wertebereich von 0 bis 2pi
        #if Ψ[:,i] > np.pi:
        #    Ψ[:,i] = Ψ[:,i] - 2*np.pi #wertebereich von -pi bis pi, funktioniert nicht mit arrays
        η[:,i] = η[:,i-1] + 1/6 * (k1+ 2*k2 + 2*k3 + k4)
    gain[index] = np.mean(η[:,-1]-η[:,1])

    #Plotten von a und b
    if k == 0 or k == 0.003:
        plt.plot(Ψ, η, ",", color="blue")
        plt.title(r"$Bewegung der Elektronen im Phasenraum, η_0={}$".format(k))
        plt.xlabel(r"$ Ψ/ rad$")
        plt.ylabel(r"$η$")
        plt.savefig(r"build/Blatt8_Aufgabe3_η{}_Feline,Hannah.pdf".format(k))
        plt.clf()

#Plotten der c
plt.plot(np.linspace(-0.01, 0.01, 21), gain)
plt.title(r"Gainkurve")
plt.xlabel(r" Energieänderung")
plt.ylabel(r"$η_0$")
plt.savefig(r"build/Blatt8_Aufgabe3c_Feline,Hannah.pdf")

