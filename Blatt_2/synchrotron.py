import numpy as np
import matplotlib.pyplot as plt

# define constants
SPEED_OF_LIGHT = 299792458 # meter per second
CIRCUMFERENCE = 100 # meter
ELECTRON_SPEED = 0.99*SPEED_OF_LIGHT # meter per second

# compute time constants
t_pass = CIRCUMFERENCE/ELECTRON_SPEED # time for one iteration
t_stop = 2*t_pass # end of simulation after two iterations
electron_radius = CIRCUMFERENCE/(2*np.pi) # radius of electron orbit
electron_frequency = 2*np.pi/t_pass # electron angular frequency
t_emission = t_pass/72 # emit light every 5 degrees, 4.68ns

def simulate_synchrotron(time):
    '''Simulate synchrotron emission by one circling electron.
    INPUT   circling time of the electron in seconds'''
    check_input(time)  
    fig, ax = plt.subplots(figsize=[6.4, 6.4]) # sizes in inches, default [6.4, 4.8]

    # draw electron motion: center in origin, specify color and transparency
    electron_curvature = plt.Circle((0, 0), electron_radius,
        linewidth=2 , edgecolor='C2', facecolor='None', alpha=0.7)
    ax.add_artist(electron_curvature) # draws circle 'electron_curvature' to axes 'ax'

    # draw electron
    x = electron_radius*np.cos(electron_frequency*time)
    y = electron_radius*np.sin(electron_frequency*time)
    ax.plot(x, y, 'C2o')

    # compute synchrotron emission
    no_of_emissions = time/t_emission
    # np.arange creates an array of evenly spaced elements within a given interval
    # here: start=0, stop=time, step=t_emission
    radiation_time = np.arange(0, time, t_emission) # every few ns a new radiation emission occurs
    # position of emission center
    radiation_x = electron_radius*np.cos(electron_frequency*radiation_time)
    radiation_y = electron_radius*np.sin(electron_frequency*radiation_time)
    # emission circle radius
    radiation_radius = SPEED_OF_LIGHT*(time-radiation_time)

    # draw synchrotron emission
    for i in range(len(radiation_time)):
        emission_shape = plt.Circle((radiation_x[i], radiation_y[i]), radiation_radius[i],
            linewidth=0.1, edgecolor='C0', facecolor='None')
        ax.add_artist(emission_shape)

    # styles
    ax.set(xlim=(-105, 105), ylim=(-105, 105), xlabel="x", ylabel="y",
        title="Umlaufzeit {:6.1f} ns".format(time*1e9)) # inserts time parameter in title
    # title format: 6 digits in total (including decimal point). One digit after decimal point
    # f represents number as float and not e.g. integer
    
    # help lines
    ax.axhline(y=0, linestyle='--', color='black', linewidth=0.5)
    ax.axvline(x=0, linestyle='--', color='black', linewidth=0.5)

    # save output
    fig.savefig('synchrotron-{:06.1f}.pdf'.format(time*1e9))
    fig.clf()


def check_input(time):
    if time > 10*t_pass:
        raise ValueError("Time exceeds recommended simulation time. Too high values may cause great computational time.")


simulate_synchrotron(0.5*t_pass)
simulate_synchrotron(1*t_pass)
simulate_synchrotron(2*t_pass)