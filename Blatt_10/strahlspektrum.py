import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as const
from scipy import signal 
from scipy.fftpack import fft

#define constants
T0 = 384*10**-9 #ns Umlaufzeit bzw. zeitlicher Abstand, zu dem das Einzelelektronenpaket im Soeicherring 
                #die punktförmige Antenne passiert
N = 1000        #Anzahl Umläufe im Elektronenspeicherring
τ = 0.5*10**-9 #Standardabweichung 
n = np.linspace(0,384000,num = 1000)
t = np.linspace(0.5*10**-9,5*10**-9, num = 1000)

#######################################################################
#######################################################################
#task a) 

#compute constants
f0 = 1/T0       #1/s Umlauffrequenz
# --> f0 = 2604166.66 1/s
ω0 = 2*np.pi/T0
# --> ω0 = 16362461 1/s
U = const.c*T0  #m Umfang des Elektronenspeicherrings
# --> U = 115.12 m

print('Die Werte U=115.12 m und f=2604166.66 1/s sprechen für das DELTA an der TU Dortmmund')

#######################################################################
#######################################################################

#task b)

#Signal pro Umlauf
def signal_an_antenne():
    ''' normalverteiltes Signal mit Standardabweichung τ, dass an Antenne alle 384ns ankommt'''
    gaussian_numbers = np.random.normal(0, τ, size=3000)
    bins = np.arange(0,384) + 0.5 - 384/2 #array([-191.5, ..., 191.5])
    histogramm = np.histogram(gaussian_numbers, bins=bins)
    data_y = np.array(histogramm[0], dtype=float)
    data_y = data_y/np.max(data_y)
    return data_y
print(signal_an_antenne(), signal_an_antenne().shape)

bins_plot = np.arange(0,383) + 1 - 384/2 #array([-191, ...., 191])
fig1,ax = plt.subplots()
ax.plot(bins_plot, signal_an_antenne())
ax.set_title(r"Signal nach einer Periode")
ax.set_xlabel(r"$t$ in $ns$")
ax.set_ylabel(r"Signal")
fig1.savefig(r"gaussfunktion.pdf")

#gesamtes Signal über 1000 Umläufe hinweg
a = signal_an_antenne()
signal_end = np.tile(a,1000)
print(signal_end, signal_end.shape)
#time = np.arange(1,383000+1) - 384/2
time = np.arange(0, len(signal_end), dtype =float)
print(time)

fig2,ax = plt.subplots()
ax.plot(time,signal_end)
ax.set_title(r"Signal nach 1000 Perioden")
ax.set_xlabel(r"$t$ in $ns$")
ax.set_ylabel(r"Signal")
fig2.savefig(r"Signal_Zeit.pdf")

#FFT
signal_fft = np.fft.fft(signal_end[0:191500])
freq = np.fft.fftfreq(int(time.shape[-1]/2))

fig3,ax = plt.subplots()
ax.plot(freq, signal_fft.real)
ax.set_title(r"FFT-Signal nach 1000 Perioden")
ax.set_xlabel(r"ω in $1/ns$")
ax.set_ylabel(r"FFT-Signal")
fig3.savefig(r"Signal_Zeit_FFT.pdf")