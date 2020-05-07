import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as const 


a = 1
t = 1 #delta t
c = const.c
beta = np.array([0.01,0])[:, np.newaxis] #anders als in der Aufgabenstellung
alpha = np.array(np.linspace(0, np.pi * 2, 360))
B = np.array([0,0]) #Punkt P in den Ursprung legen

d = np.cos(alpha) /a
b = np.sin(alpha) /a
r = c * t
P_ret = np.array([d + beta[0]*r,b])
n = (P_ret.T - B /r).T

E = (1-beta[0]**2)*(n-beta)/(r**2*(1-np.dot(beta.T, n))**3)
print(E)

