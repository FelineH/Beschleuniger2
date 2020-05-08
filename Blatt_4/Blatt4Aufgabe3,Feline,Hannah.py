import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 

#definiere die konstanten und startwerte
a = 1 #parameter aus der skizze
t = 1 #delta t aus der skizze
c = const.c #lichtgeschwindigkeit in m/s
beta = np.array([0.9,0])[:, np.newaxis] #beta wird erzeugt aber anders als in der Aufgabenstellung, np.newaxis ist für später damit man es mit einer matrix multiplizieren kann
alpha = np.array(np.linspace(0, np.pi * 2, 360)) #array mit allen winkeln wird erstellt
B = np.array([0,0])[:, np.newaxis] #Punkt B in den Ursprung legen

#längen und vektoren aus der skizze berechnen
d = np.cos(alpha) *a 
b = np.sin(alpha) *a
r = c*t #hier passt die größenordung noch nicht da alle anderen längen in willkürlichen einheiten sind
P_ret = np.array([d + beta[0]*r,b]) #P_ret entspricht P'
n = ((P_ret - B) /r) 
n = n/np.linalg.norm(n) #normieren
r = r/np.linalg.norm(n) #r mit dem gleichen faktor normieren, damit die länge in den gleichen einheiten ist wie der rest

#E-Feld berechen, Transponieren damit die größen zusammen passen
E = (1-beta[0]**2)*(n-beta)/(r**2*(1-np.dot(beta.T,n))**3)

fig, ax = plt.subplots()
plt.quiver(0,0,E[0,:],E[1,:], width=0.0008,)
ax.set_title('Elektrisches Feld')
ax.set_xlabel(r"$x $")
ax.set_ylabel(r"$y $")
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
fig.savefig('feld.pdf')
fig.clf()
print(E.shape)
