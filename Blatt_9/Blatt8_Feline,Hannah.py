import numpy as np
import matplotlib.pyplot as plt

#Konstanten
λ_u = 0.027 #m
K = 1.19
ɣ = 2000
m_e = 9.1e-31 #kg
c = 3e8 #m/s
e = 1.6e-19 #C
μ_0 = 1.257e-6 #N/A²
n_e = 1e22 #m⁻³

N = 1000 #anzahl elektronen
x = 0.001 #m, wegschritte

a1 = 4 * np.pi / λ_u
a2 = - e * K / (2 * m_e * c**2 * ɣ**2)
a3 = e * μ_0 * c**2 * K * n_e / (2 * ɣ)


leistung = np.zeros(6000) #für c


#Differentialgleichungen
def dΨ(η, a1,):
    return a1 * η

def dη(E, Ψ, a2):
    return a2 * np.real(E * np.exp(1j * Ψ)) 

def dE(Ψ, a3,):
    return a3 * 1/N *np.sum(np.exp(-1j * Ψ))


#Startwerte definieren
Ψ = np.zeros((N, 7))
Ψ[:,0] = np.random.uniform(0, 4*np.pi, N) #gleichverteilte Phase, von 0 bis 2pi, weil es schöner dargestellt wird, ändert nichts an berechnung
Ψakt = Ψ[:,0]
η = np.zeros((N, 7))
η[:,0] = np.random.normal(scale= 0.001, size=N)
ηakt = η[:,0]
E = np.zeros(7, dtype=np.complex)
E[0] = np.complex(1e6, 0j)
Eakt = E[0]


#Runge-Kutta-Verfahren
for n in range(1, 7):
    for i in range(1, 1001):

        j1 = x * dΨ(ηakt, a1,)
        k1 = x * dη(Eakt,Ψakt, a2)
        m1 = x * dE(Ψakt, a3)
        j2 = x * dΨ(ηakt+k1/2, a1,)
        k2 = x * dη(Eakt+m1/2,Ψakt+j1/2, a2)
        m2 = x * dE(Ψakt+j1/2, a3)
        j3 = x * dΨ(ηakt+k2/2, a1,)
        k3 = x * dη(Eakt+m2/2,Ψakt+j2/2, a2)
        m3 = x * dE(Ψakt+j2/2, a3)
        j4 = x * dΨ(ηakt+k3, a1,)
        k4 = x * dη(Eakt+m3,Ψakt+j3, a2)
        m4 = x * dE(Ψakt+j3, a3)

        Ψakt = Ψakt + 1/6 * (j1+ 2*j2 + 2*j3 + j4)
        Ψakt = Ψakt%(4*np.pi) #Wertebereich von 0 bis 2pi
        #if Ψ[:,i] > np.pi:
        #    Ψ[:,i] = Ψ[:,i] - 2*np.pi #wertebereich von -pi bis pi, funktioniert nicht mit arrays
        ηakt = ηakt + 1/6 * (k1+ 2*k2 + 2*k3 + k4)
        Eakt = Eakt + 1/6 * (m1+ 2*m2 + 2*m3 + m4)
        leistung[1000*(n-1)+(i-1)] = np.real(Eakt)**2 + np.imag(Eakt)**2

    Ψ[:,n] = Ψakt
    η[:,n] = ηakt
    E[n] = Eakt

    plt.plot(Ψakt, ηakt, ".", color="blue")
    plt.title(r"$Bewegung der Elektronen im Phasenraum, s={}m$".format(n))
    plt.xlabel(r"$ Ψ/ rad$")
    plt.ylabel(r"$η$")
    plt.savefig(r"Blatt8_Aufgabe3_s{}_Feline,Hannah.pdf".format(n))
    plt.clf()

plt.plot(np.linspace(0,6, num=6000), leistung)
plt.title(r"Leistung, E_0={},n_e={}".format(E[0],n_e))
plt.xlabel(r"$ s/ m")
plt.ylabel(r"$η$")
plt.savefig(r"Leistung,E_0={},n_e={}.pdf".format(E[0],n_e))






