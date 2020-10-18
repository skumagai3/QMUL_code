# Plotting RHS of eq. 15 of Perico et al. 2019
# (Updated equations)
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.offsetbox import AnchoredText

y = np.linspace(0.00000001,5,8192)

# User input
choice = input('Symmetron or f(R)? Type S or f. ')
choice = str(choice)
if choice == 'f':
    exponent = float(input('Choose |f_R0| value (10**-4, -5, or -6): '))
    f_R0_choice = 10**(-exponent)
    print(f_R0_choice)
if choice == 'S':
    z_ssb_choice = int(input('Choose z_ssb value (0 to 3): '))
    print(z_ssb_choice)

sigma_0 = 0.811
delta_0 = -0.8
alpha_0 = -1.52
beta_0 = 3.0
c_0 = 3.0
c_1 = 0.3
G_0 = 0.3


def gamma_fR(f_R0,gamma_1=2.0,gamma_2=0.1,gamma_3=1e-5):
    return gamma_1+gamma_2*np.log10(gamma_3+abs(f_R0))

def gamma_sym(z_SSB,gamma_4=2.2,gamma_5=0.15,gamma_6=0.15):
    return gamma_4+gamma_5*(z_SSB/(1+gamma_6*math.exp(-(z_SSB**2))))

# VOID DENSITY PROFILE
def RHS_fR(f_R0,sigma_0=sigma_0,delta_0=delta_0,alpha_0=alpha_0,
           beta_0=beta_0,c_0=c_0,c_1=c_1,G_0=G_0):
    alpha = gamma_fR(f_R0)-alpha_0*sigma_0
    beta = beta_0*alpha
    c = (c_0*sigma_0**(-c_1))**(1/beta)
    G = gamma_fR(f_R0)+G_0*np.log10(sigma_0)
    RHS = delta_0*((1-G*(y*c)**alpha)/(1+(y*c)**beta))
    return RHS

def RHS_SYM(z_SSB,sigma_0=sigma_0,delta_0=delta_0,alpha_0=alpha_0,
           beta_0=beta_0,c_0=c_0,c_1=c_1,G_0=G_0):
    alpha = gamma_sym(z_SSB)-alpha_0*sigma_0
    beta = beta_0*alpha
    c = (c_0*sigma_0**(-c_1))**(1/beta)
    G = gamma_sym(z_SSB)+G_0*np.log10(sigma_0)
    RHS = delta_0*((1-G*(y*c)**alpha)/(1+(y*c)**beta))
    return RHS

# Plotting
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlabel(r'Radius')
ax.set_ylabel(r'$\frac{\rho_{v}(r)}{\rho_{m}}$')
ax.axhline(y=0,linestyle='dotted',c='k')
if choice == 'S':
    ax.plot(y,RHS_SYM(z_SSB=z_ssb_choice))
    ax.set_title('Density Profile of a Cosmic Void (Symmetron)')
    anchored_text = AnchoredText(r'$z_{SSB}=$'+str(z_ssb_choice)
	+'\n'+r'$\sigma=$'+str(sigma_0)
    +'\n'+r'$\delta_{0}=$'+'%s' % float('%.3g' % delta_0)
    +'\n'+r'$\alpha_{0}=$'+'%s' % float('%.3g' % alpha_0)
    +'\n'+r'$\beta_{0}=$'+'%s' % float('%.3g' % beta_0)
    +'\n'+r'$c_{0}=$'+'%s' % float('%.3g' % c_0)
    +'\n'+r'$c_{1}=$'+'%s' % float('%.3g' % c_1)
    +'\n'+r'$G_{0}=$'+'%s' % float('%.3g' % G_0),loc=4)
else:
    ax.plot(y,RHS_fR(f_R0=f_R0_choice))
    ax.set_title('Density Profile of a Cosmic Void (f(R))')
    anchored_text = AnchoredText(r'$|f_{R0}|=$'+str(f_R0_choice)
	+'\n'+r'$\sigma=$'+str(sigma_0)
    +'\n'+r'$\delta_{0}=$'+'%s' % float('%.3g' % delta_0)
    +'\n'+r'$\alpha_{0}=$'+'%s' % float('%.3g' % alpha_0)
    +'\n'+r'$\beta_{0}=$'+'%s' % float('%.3g' % beta_0)
    +'\n'+r'$c_{0}=$'+'%s' % float('%.3g' % c_0)
    +'\n'+r'$c_{1}=$'+'%s' % float('%.3g' % c_1)
    +'\n'+r'$G_{0}=$'+'%s' % float('%.3g' % G_0),loc=4)
ax.add_artist(anchored_text)
if choice == 'S':
    fig.savefig('new_profile_sym_z'+str(z_ssb_choice)+'.png')
    print('Figure has been saved.')
if choice == 'f':
    fig.savefig('new_profile_f(R)_f'+str(f_R0_choice)+'.png')
    print('Figure has been saved.')
plt.show()