# Plotting RHS of eq. 15 of Perico et al. 2019
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

sym = np.linspace(1,3,2048)
fr = np.linspace(10**-6,10**-4,2048)

# Choose constants
alpha0=3.0
alpha1=0.3
alpha2=10**-5
alpha3 = 0.3
beta0 = 7.0
c0 = 0.3
c1 = 0.3
c2 = 0.3
G0 = 0.3
G1 = 0.3
G2 = 10**-5
G3 = 0.3
G4 = 0.3

# User input
choice = input('Symmetron or f(R)? Type S or f. ')
choice = str(choice)

# NEW Sigma value:
sigma0 = 0.811
# Gamma functional form for Symmetron
def gamma_sym(gamma_0,gamma_1,gamma_2,z_ssb):
    return gamma_0 + gamma_1*(z_ssb/(1+gamma_2*math.exp(-z_ssb**2)))
# Gamma functional form for f(R)
def gamma_R(gamma_0,gamma_1,gamma_2,f_R0):
	return gamma_0 + gamma_1*np.log10(gamma_2+f_R0)

# alpha fxns
def alpha_sym(alpha_0=alpha0,alpha_1=alpha1,alpha_2=alpha2,alpha_3=alpha3,zssb=2):
    return gamma_sym(alpha_0,alpha_1,alpha_2,zssb)+alpha_3*sigma0
def alpha_R(alpha_0=alpha0,alpha_1=alpha1,alpha_2=alpha2,alpha_3=alpha3,f_R0=10**-5):
    return gamma_R(alpha_0,alpha_1,alpha_2,f_R0)+alpha_3*sigma0
# beta fxns
def beta_sym(beta_0=beta0,zssb=2):
    return beta_0*alpha_sym(zssb=zssb)
def beta_R(beta_0=beta0,f_R0=10**-5):
    return beta_0*alpha_R(f_R0=f_R0)
# c fxns
def c_sym(c_0=c0,c_1=c1,c_2=c2,zssb=2):
    return (c_0+c_1*sigma0**c_2)**(-1/beta_sym(zssb=zssb))
def c_R(c_0=c0,c_1=c1,c_2=c2,f_R0=10**-5):
    return (c_0+c_1*sigma0**c_2)**(-1/beta_R(f_R0=f_R0))
# G fxns
def G_sym(G_0=G0,G_1=G1,G_2=G2,G_3=G3,G_4=G4,zssb=2):
    return gamma_sym(G_0,G_1,G_2,zssb)+G_3*sigma0+G_4*sigma0**2
def G_R(G_0=G0,G_1=G1,G_2=G2,G_3=G3,G_4=G4,f_R0=10**-5):
    return gamma_R(G_0,G_1,G_2,f_R0)+G_3*sigma0+G_4*sigma0**2

# Make figure area with space for chosen values at bottom
fig, ax = plt.subplots(figsize=(8, 6))
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.3)
plt.gcf().text(0.1, 0.13,r'$\alpha_{0}=$'+str(alpha0)+', '
    +r'$\alpha_{1}=$'+str(alpha1)+', '
    +r'$\alpha_{2}=$'+str(alpha2)+', '
    +r'$\alpha_{3}=$'+str(alpha3)+', '
    +r'$\beta_{0}=$'+str(beta0)+', '
    +r'$c_{0}=$'+str(c0)+', '
    +r'$c_{1}=$'+str(c1)+', '
    +r'$c_{2}=$'+str(c2)+', '
    +'\n'+r'$G_{0}=$'+str(G0)+', '
    +r'$G_{1}=$'+str(G1)+', '
    +r'$G_{2}=$'+str(G2)+', '
    +r'$G_{3}=$'+str(G3)+', '
    +r'$G_{4}=$'+str(G4))

# Making graphs
if choice == 'S':
    alpha_values = []
    beta_values = []
    c_values = []
    G_values = []
    for value in sym:
        alpha_values.append(alpha_sym(zssb=value))
        beta_values.append(beta_sym(zssb=value))
        c_values.append(c_sym(zssb=value))
        G_values.append(G_sym(zssb=value))
    ax.set_xlabel(r'$z_{SSB}$')
    ax.set_ylabel(r'Function Value')
    ax.set_title(r'$Functions\ vs\ Parameter\ x=z_{SSB}$')
    ax.plot(sym,alpha_values,label=r'$\alpha$')
    ax.plot(sym,beta_values,label=r'$\beta$')
    ax.plot(sym,c_values,label='c')
    ax.plot(sym,G_values,label='G')
if choice == 'f':
    alpha_values = []
    beta_values = []
    c_values = []
    G_values = []
    for value in fr:
        alpha_values.append(alpha_R(f_R0=value))
        beta_values.append(beta_R(f_R0=value))
        c_values.append(c_R(f_R0=value))
        G_values.append(G_R(f_R0=value))
    ax.set_xlabel(r'$|f_{R0}|$')
    ax.set_ylabel(r'Function Value')
    ax.set_title(r'$Functions\ vs\ Parameter\ x=|f_{R0}|$')
    ax.plot(fr,alpha_values,label=r'$\alpha$')
    ax.plot(fr,beta_values,label=r'$\beta$')
    ax.plot(fr,c_values,label='c')
    ax.plot(fr,G_values,label='G')
plt.legend()
# Automatically save figure
if choice == 'S':
    fig.savefig('params_sym.png')
    print('Figure has been saved.')
if choice == 'f':
    fig.savefig('params_fR.png')
    print('Figure has been saved.')

plt.show()