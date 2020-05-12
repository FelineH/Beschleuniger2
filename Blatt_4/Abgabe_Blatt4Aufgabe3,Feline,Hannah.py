import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

#define constants
 
#a = 1 #parameter aus der skizze
#c = const.c #lichtgeschwindigkeit in m/s

beta = np.array([0.01, 0.50, 0.99]) #konstante Geschwindigkeit des Teilchen in z-Richtung, wobei beta 0.01, 0.5 und 0.99 sein soll 
alpha = np.linspace(0, 2*np.pi, 360) #Array mit allen Winkeln von 0 bis 2pi in 360 Schritten

length = 10 #meter #length of the plotted arrows



#define functions

def viewer(alpha): #(meter, meter)^T als 2D Spaltenvektor 
    """ Beobachter im Punkt B, der sich um Punkt P herumbewegt
    """
    B = np.array([np.cos(alpha), np.sin(alpha)])
    return B
    #print (type(B), B.shape) <class 'numpy.ndarray'> (2, 360) 2 Zeilen mit 360 Spalten

#Wende p-q-Formel an, um r zu errechnen
def r(d, b):
    """ Berechne Skalar r mithilfe der p-q-Formel und der Relationen aus der Skizze 
    mit r**2 = b**2 + (d+beta*r)**2, die nach r**2 + (2*beta*d/((beta**2)-1))*r + (b**2+d**2)/((beta**2)-1) umgestellt wird 
    """
    p = 2*d*beta[2]/(beta[2]**2-1) # p entspricht 
    q = (b**2+d**2)/(beta[1]**2-1)
    r = -(p/2) + np.sqrt((p/2)**2 - q)
    return r
    #print(type(r), r.shape, r) <class 'numpy.ndarray'> (360,) mit 360 Spalten und 1 Zeile, heißt 1d Array bzw. Zeilenvektor

def n_vec(B):
    """ Berechne normierten Einheitsvektor n_norm von P' in Richtung B 
    """
    P_ret = np.array([-beta[2] * r(B[0], B[1]), np.zeros(360)]) #P_ret entspricht P' in der Skizze, wobei P# dem retadierten Ort der Ladung entspricht 
    n = (B - P_ret) 
    #print(n)
    n_norm = n /np.linalg.norm(n) #normieren
    return n_norm
    #return P_ret
print(n_vec(viewer(alpha)), n_vec(viewer(alpha)).shape) #(360,2) mit 360 Spalten und 2 Zeilen (Dimension entspricht P_ret)
    
#Richtung des E-Vektors berechnen in Abhängigkeit des Einheitsvektors n_vec
direction = length*n_vec(viewer(alpha))

print('Richtung', direction, direction.shape) #(360,2)
print('beobachter(alpha)[0]', viewer(alpha)[0], viewer(alpha)[0].shape) #(360,)
print('richtung[0]', direction[0], direction[0].shape) #(360,)




#Plot arrows

fig, ax = plt.subplots()
x = np.any(viewer(alpha)[0])
y = np.any(viewer(alpha)[1])
u = np.any(direction[0])
v = np.any(direction[1])
ax.quiver(x, y, u, v, linewidth=0.1)
##first two parameters (here: np.any(viewer(alpha)[0]), np.any(viewer(alpha)[1])) define the arrow location, 
#the third and fourth parameters defíne the arrow direction (here: the direction of the E-vector)
#the arrays can be 1D or 2D like 
ax.axis('equal')
ax.plot([0],[0],'ro') #Ursprungspunkt, von dem aus im Abstand a=1 in 1° Grad Schritten des Winkels alpha E-Vektoren gezeichnet werden sollen
ax.set_title('Elektrisches Feld')
ax.set_xlabel(r"$x $")
ax.set_ylabel(r"$y $")
ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
fig.savefig('build/Feld_beta_0.99.pdf')
fig.clf()


