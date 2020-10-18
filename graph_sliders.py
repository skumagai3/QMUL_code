# Plotting RHS of eq. 15 of Perico et al. 2019
import math
import pylab
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.constants import c
from matplotlib.offsetbox import AnchoredText
from matplotlib.widgets import Slider, Button, RadioButtons
# Create radius axis
y = np.linspace(0,3,2048)

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

# EQUATION
def RHS_SYM(z_ssb=2,delta0=-0.8,alpha3=0.3,beta0=7,c0=0.3,c1=0.3,c2=0.3,G3=0.3,G4=0.3):
    alpha = gamma_sym(3,0.3,10**-5,z_ssb)+alpha3*sigma0
    beta = beta0*alpha
    c = (c0+c1*sigma0**c2)**(-1/beta)
    G = gamma_sym(0.3,0.3,10**-5,z_ssb)+G3*sigma0+G4*sigma0**2
    RHS = delta0*(1-G*(y*c)**alpha)/(1+(y*c)**beta)
    return RHS
def RHS_FR(f_R=-5,delta0=-0.8,alpha3=0.3,beta0=7,c0=0.3,c1=0.3,c2=0.3,G3=0.3,G4=0.3):
    alpha = gamma_R(3,0.1,10**-5,10**f_R)+alpha3*sigma0
    beta = beta0*alpha
    c = (c0+c1*sigma0**c2)**(-1/beta)
    G = gamma_R(0.3,0.1,10**-5,10**f_R)+G3*sigma0+G4*sigma0**2
    RHS = delta0*(1-G*(y*c)**alpha)/(1+(y*c)**beta)
    return RHS

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
ax.set_xlabel(r'Radius')
ax.set_ylabel(r'$\frac{\rho_{v}(r)}{\rho_{m}}$')
ax.axhline(y=0,linestyle='dotted',c='k')
fig.subplots_adjust(top=0.98)
fig.subplots_adjust(bottom=0.6)

# Choice between Symmetron and f(R)
if choice == 'S':
    [line] = ax.plot(y, RHS_SYM(), linewidth=2, color='blue')
else:
    [line] = ax.plot(y, RHS_FR(), linewidth=2, color='blue')

# Define axes area and make sliders for either MG model
if choice == 'S':
    S_slider_ax = fig.add_axes([0.25, 0.5, 0.65, 0.03])
    S_slider = Slider(S_slider_ax, r'$z_{ssb}$',1,3,valinit=2,valstep=1.0)
else:
    FR_slider_ax = fig.add_axes([0.25, 0.5, 0.65, 0.03])
    FR_slider = Slider(FR_slider_ax, r'$f(R)$',-6,-4,valinit=-5,valstep=1.0)

# Make sliders for delta,alpha,cs,Gs
delta0_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
delta0_slider = Slider(delta0_slider_ax, r'$\delta_{0}$', -1.0, -0.1, valinit=-0.8)

alpha3_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
alpha3_slider = Slider(alpha3_slider_ax, r'$\alpha_{3}$', 0.01, 2, valinit=0.3)

beta0_slider_ax = fig.add_axes([0.25, 0.2, 0.65, 0.03])
beta0_slider = Slider(beta0_slider_ax, r'$\beta_{0}$', 1, 10, valinit=7)

c0_slider_ax = fig.add_axes([0.25, 0.25, 0.65, 0.03])
c0_slider = Slider(c0_slider_ax, r'$c_{0}$', 0.01, 2, valinit=0.3)
c1_slider_ax = fig.add_axes([0.25, 0.3, 0.65, 0.03])
c1_slider = Slider(c1_slider_ax, r'$c_{1}$', 0.01, 2, valinit=0.3)
c2_slider_ax = fig.add_axes([0.25, 0.35, 0.65, 0.03])
c2_slider = Slider(c2_slider_ax, r'$c_{2}$', 0.01, 2, valinit=0.3)

G3_slider_ax = fig.add_axes([0.25, 0.4, 0.65, 0.03])
G3_slider = Slider(G3_slider_ax, r'$G_{3}$', 0.01, 2, valinit=0.3)
G4_slider_ax = fig.add_axes([0.25, 0.45, 0.65, 0.03])
G4_slider = Slider(G4_slider_ax, r'$G_{4}$', 0.01, 2, valinit=0.3)

# Define an action for modifying the line when any slider's value changes
def sliders_on_changed_SYM(val):
    line.set_ydata(RHS_SYM(z_ssb=S_slider.val,delta0=delta0_slider.val,alpha3=alpha3_slider.val,
    beta0=beta0_slider.val,c0=c0_slider.val,c1=c1_slider.val,c2=c2_slider.val,G3=G3_slider.val,
    G4=G4_slider.val))
    fig.canvas.draw_idle()
def sliders_on_changed_FR(val):
    line.set_ydata(RHS_FR(f_R=FR_slider.val,delta0=delta0_slider.val,alpha3=alpha3_slider.val,
    beta0=beta0_slider.val,c0=c0_slider.val,c1=c1_slider.val,c2=c2_slider.val,G3=G3_slider.val,
    G4=G4_slider.val))
    fig.canvas.draw_idle()

if choice == 'S':
    S_slider.on_changed(sliders_on_changed_SYM)
    alpha3_slider.on_changed(sliders_on_changed_SYM)
    delta0_slider.on_changed(sliders_on_changed_SYM)
    beta0_slider.on_changed(sliders_on_changed_SYM)
    c0_slider.on_changed(sliders_on_changed_SYM)
    c1_slider.on_changed(sliders_on_changed_SYM)
    c2_slider.on_changed(sliders_on_changed_SYM)
    G3_slider.on_changed(sliders_on_changed_SYM)
    G4_slider.on_changed(sliders_on_changed_SYM)
else:
    FR_slider.on_changed(sliders_on_changed_FR)
    alpha3_slider.on_changed(sliders_on_changed_FR)
    delta0_slider.on_changed(sliders_on_changed_FR)
    beta0_slider.on_changed(sliders_on_changed_FR)
    c0_slider.on_changed(sliders_on_changed_FR)
    c1_slider.on_changed(sliders_on_changed_FR)
    c2_slider.on_changed(sliders_on_changed_FR)
    G3_slider.on_changed(sliders_on_changed_FR)
    G4_slider.on_changed(sliders_on_changed_FR)

# Add a button for resetting the parameters
axis_color = 'lightgoldenrodyellow'
reset_button_ax = fig.add_axes([0.1, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    alpha3_slider.reset()
    delta0_slider.reset()
    beta0_slider.reset()
    c0_slider.reset()
    c1_slider.reset()
    c2_slider.reset()
    G3_slider.reset()
    G4_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)
# Create save button that saves density profile NOT sliders
save_button_ax = fig.add_axes([0.1, 0.1, 0.1, 0.04])
save_button = Button(save_button_ax,'Save',color=axis_color,hovercolor='0.9')
def save_button_on_clicked(mouse_event):
    print('Figure has been saved.')
    fig1 = plt.figure(2,figsize=(7,5))
    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel(r'Radius')
    ax1.set_ylabel(r'$\frac{\rho_{v}(r)}{\rho_{m}}$')
    ax1.axhline(y=0,linestyle='dotted',c='k')
    if choice == 'S':
        ax1.plot(y, RHS_SYM(z_ssb=S_slider.val,delta0=delta0_slider.val,
                alpha3=alpha3_slider.val,beta0=beta0_slider.val,c0=c0_slider.val,
                c1=c1_slider.val,c2=c2_slider.val,G3=G3_slider.val,
                G4=G4_slider.val), linewidth=2, color='blue')
        anchored_text = AnchoredText(r'$z_{SSB}=$'+str(S_slider.val)
        +'\n'+r'$\sigma=$'+str(sigma0)
        +'\n'+r'$\delta_{0}=$'+'%s' % float('%.3g' % delta0_slider.val)
        +'\n'+r'$\alpha_{3}=$'+'%s' % float('%.3g' % alpha3_slider.val)
        +'\n'+r'$\beta_{0}=$'+'%s' % float('%.3g' % beta0_slider.val)
        +'\n'+r'$c_{0}=$'+'%s' % float('%.3g' % c0_slider.val)
        +'\n'+r'$c_{1}=$'+'%s' % float('%.3g' % c1_slider.val)
        +'\n'+r'$c_{2}=$'+'%s' % float('%.3g' % c2_slider.val)
        +'\n'+r'$G_{3}=$'+'%s' % float('%.3g' % G3_slider.val)
        +'\n'+r'$G_{4}=$'+'%s' % float('%.3g' % G4_slider.val), loc=4)
        ax1.add_artist(anchored_text)
        ax1.set_title('Density Profile of a Cosmic Void (Symmetron)')
        fig1.savefig('profile_sym_zssb='+str(S_slider.val)+'.png')
    else:
        ax1.plot(y, RHS_FR(f_R=FR_slider.val,delta0=delta0_slider.val,
                alpha3=alpha3_slider.val,beta0=beta0_slider.val,c0=c0_slider.val,
                c1=c1_slider.val,c2=c2_slider.val,G3=G3_slider.val,
                G4=G4_slider.val), linewidth=2, color='blue')
        anchored_text = AnchoredText(r'$|f_{R0}|=$'+str(FR_slider.val)
        +'\n'+r'$\sigma=$'+str(sigma0)
        +'\n'+r'$\delta_{0}=$'+'%s' % float('%.3g' % delta0_slider.val)
        +'\n'+r'$\alpha_{3}=$'+'%s' % float('%.3g' % alpha3_slider.val)
        +'\n'+r'$\beta_{0}=$'+'%s' % float('%.3g' % beta0_slider.val)
        +'\n'+r'$c_{0}=$'+'%s' % float('%.3g' % c0_slider.val)
        +'\n'+r'$c_{1}=$'+'%s' % float('%.3g' % c1_slider.val)
        +'\n'+r'$c_{2}=$'+'%s' % float('%.3g' % c2_slider.val)
        +'\n'+r'$G_{3}=$'+'%s' % float('%.3g' % G3_slider.val)
        +'\n'+r'$G_{4}=$'+'%s' % float('%.3g' % G4_slider.val), loc=4)
        ax1.add_artist(anchored_text)
        ax1.set_title('Density Profile of a Cosmic Void (f(R))')
        fig1.savefig('profile_f(R)_1e'+str(FR_slider.val)+'.png')
    plt.show()
save_button.on_clicked(save_button_on_clicked)
plt.show()