import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c

import os 
if not os.path.isdir('build'):
    os.mkdir('build')

a_begin = 0.5
a_end = 0.6
N_a = 1
alpha_begin = -180
alpha_end = 180
N_alpha = alpha_end - alpha_begin
B_W = np.asarray([0, 1], np.float64)
plot_window = 5
factor = 1


# "retarded" position (:P)
def cal_P_prime_W(B_W, beta_W):
    d, b = B_W
    beta_W = np.linalg.norm(beta_W)
    beta = np.linalg.norm(beta_W)
    p = 2*beta**2*d/(1-beta**2)
    q = -1*(beta**2*d**2+beta**2*b**2)/(1-beta**2)
    d_prime_m = -p/2 - np.sqrt((p/2)**2-q)
    d_prime_p = -p/2 + np.sqrt((p/2)**2-q)
    return np.asarray([d_prime_m, 0], np.float64)


# coulomb term
def cal_E_W(beta_W, n_W, r):
    beta = np.linalg.norm(beta_W)
    return (1-beta**2) * (n_W - beta_W) / (r**2*(1-np.einsum('w,w->', n_W, beta_W))**3)


# calculations
analysis_results = {}
for beta in [0.01, 0.5, 0.99]:
    print('running analysis for beta = {:.2f}'.format(beta))
    B_Ws = []
    E_Ws = []
    alphas = []
    for alpha in np.linspace(alpha_begin, alpha_end, N_alpha):
        for a in np.linspace(a_begin, a_end, N_a):
            tan = np.tan(2*np.pi/360*alpha)
            b = a*np.sin(2*np.pi/360*alpha)
            d = a*np.cos(2*np.pi/360*alpha)
            beta_W = np.asarray([beta, 0], np.float64)
            B_W = np.asarray([d, b], np.float64)
            P_prime_W = cal_P_prime_W(B_W, beta_W)
            r = np.linalg.norm(P_prime_W - B_W)
            n_W = (B_W-P_prime_W) / r
            E_W = cal_E_W(beta_W, n_W, r)
            alphas.append(alpha)
            E_Ws.append(E_W)
            B_Ws.append(B_W)
    alpha_M = np.asarray(alphas)
    B_MW = np.asarray(B_Ws)
    E_MW = np.asarray(E_Ws)
    analysis_results.update({beta: {'B_MW': B_MW, 'E_MW': E_MW, 'alpha_M': alpha_M}})


# plots for a)
# ------------
for beta in [0.01, 0.5, 0.99]:
    print('plotting a) for beta = {:.2f}'.format(beta))
    B_MW = analysis_results[beta]['B_MW']
    E_MW = analysis_results[beta]['E_MW']
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111)
    ax.set_aspect('equal')
    plt.title(r'$\beta={:.2f}$'.format(beta))
    for B_W, E_W in zip(B_MW, E_MW):
        E_norm_W = factor*E_W  # /np.linalg.norm(E_W)*np.linalg.norm(E_W)**(1/2)
        plt.arrow(*B_W, *E_norm_W, length_includes_head=True)
    plt.plot([0], [0], 'o', color='r', label='charge')
    plt.plot([], [], color='0', label='coulomb field')
    plt.xlim((-plot_window, plot_window))
    plt.ylim((-plot_window, plot_window))
    plt.xlabel('$x$ (a.u.)')
    plt.ylabel('$y$ (a.u.)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('build/E_beta{}.png'.format(int(beta*100)), mode='png', dpi=300)
    plt.savefig('build/E_beta{}.pdf'.format(int(beta*100)), mode='pdf')


# plots for b)
# ------------
alpha_beta001_M = analysis_results[0.01]['alpha_M']
B_beta001_MW = analysis_results[0.01]['B_MW']
E_beta001_MW = analysis_results[0.01]['E_MW']
E_beta001_M = np.linalg.norm(E_beta001_MW, axis=1)
for beta in [0.5, 0.99]:
    print('plotting b) for beta = {:.2f}'.format(beta))
    alpha_M = analysis_results[beta]['alpha_M']
    B_MW = analysis_results[beta]['B_MW']
    E_MW = analysis_results[beta]['E_MW']
    E_M = np.linalg.norm(E_MW, axis=1)
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(111, projection='polar')
    ax.set_aspect('equal')
    plt.title(r'$E(\beta={:.2})/E(\beta=0.01)$ vs. $\alpha$'.format(beta))
    plt.plot(2*np.pi/360*alpha_M, E_M/E_beta001_M, '+')
    plt.tight_layout()
    plt.savefig('build/field_ratio_beta{}.png'.format(int(beta*100)), mode='png', dpi=300)
    # plt.savefig('field_quotient_beta{}.pdf'.format(int(beta*100)), mode='pdf')

