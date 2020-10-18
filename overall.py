"""
Adding ignored OmegaL0 version:
9/7 version 1
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True )
plt.rc('font', family='Times New Roman', weight='normal', size=14)
plt.rcParams['mathtext.fontset'] = 'dejavusans' #'Times New Roman'
plt.rcParams.update({'axes.labelsize': 'large'})
plt.rcParams.update({'axes.titlesize': 'large'})
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
# Constants
G = 6.674e-11           #m^3 kg^-1 s^-2
Mpc = 3.086e19			#km

# bkgd present time params
omega_m = 0.30
omega_l = 0.70
omega_r = 0.00
H0 = 70

# User inputs:
final_eta_value = float(input('Input final eta value: '))
w_0, w_a = input('Input w_0, w_a values: ').split()
w_0 = float(w_0)
w_a = float(w_a)
div_time = float(input('Input time of void deviation from bkgd universe: '))

#################################################################################
# LCDM FUNCTIONS
#################################################################################
# BACKGROUND
def lcdm_da_dt(a,t):
    return H0*np.sqrt(omega_m/a + omega_l*a**2)
def lcdm_H2(a):
    return H0**2 * (omega_m/a**3 + omega_l)
def lcdm_omg_M(a):
    return (omega_m*H0**2)/(lcdm_H2(a)*a**3)
def lcdm_omg_L(a):
    return (omega_l*H0**2)/(lcdm_H2(a))
# VOID
def lcdm_eta_fxn(a):
    A = final_eta_value/omega_l
    return A * lcdm_omg_L(a)
def lcdm_da_dt_v(a,t):
    eta_value = lcdm_eta_fxn(a)
    return H0 * np.sqrt(omega_m/a + omega_l*a**2 - 5.*eta_value*omega_m/3.)
def lcdm_H_v2(a):
    return H0**2 * (omega_m/a**3 + omega_l - 5*lcdm_eta_fxn(a)*omega_m/(3.*a**2))
def lcdm_omg_vM(a):
    return (omega_m*H0**2)/(lcdm_H_v2(a)*a**3)
def lcdm_omg_vL(a):
    return (omega_l*H0**2)/(lcdm_H_v2(a))
def lcdm_omg_vk(a):
    return -H0**2 * (5*lcdm_eta_fxn(a)*omega_m)/(3.*lcdm_H_v2(a)*a**2)


#################################################################################
# OLD CPL PARAMETRIZATION FUNCTIONS
#################################################################################
# BACKGROUND
def da_dt(a,t,w_0,w_a):
    return H0*(omega_m/a + omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a)))**(1./2.)
def H2(a,w_0,w_a):
    return H0**2*(omega_m/a**3 + omega_l*np.exp(3.*w_a*(a-1.))/a**(3.*(1.+w_0+w_a)))
def omg_M(a,w_0,w_a):
    return (omega_m*H0**2)/(a**3*H2(a,w_0,w_a))
def omg_L(a,w_0,w_a):
    return (omega_l*H0**2)/(H2(a,w_0,w_a))*(np.exp(3.*w_a*(a-1.)))/(a**(3*(1+w_0+w_a)))
# VOID
def eta_fxn(a,w_0,w_a):
    return final_eta_value * (H0**2)/(H2(a,w_0,w_a)) * (np.exp(3*w_a*(a-1.)))/(a**(3.*(1.+w_0+w_a)))
def da_dt_v(a,t,w_0,w_a):
    eta_value = eta_fxn(a,w_0,w_a)
    deriv = H0*(omega_m/a+omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a))
    -eta_value/3.*(5.*omega_m+2.*omega_l*np.exp(3.*w_a*(a-1.))/a**(3.*(w_0+w_a))))**(1./2.)
    return deriv
def H_v2(a,w_0,w_a):
    val = H0**2*(omega_m/a**3 
    + omega_l*np.exp(3.*w_a*(a-1.))/a**(3*(1+w_0+w_a)) 
    - eta_fxn(a,w_0,w_a)/(3.*a**2)*(5.*omega_m+2.*omega_l*np.exp(3.*w_a*(a-1.))/a**(3*(w_0+w_a))))
    return val
def omg_vM(a,w_0,w_a):
    return (omega_m*H0**2)/(H_v2(a,w_0,w_a)*a**3)
def omg_vL(a,w_0,w_a):
    return (omega_l*H0**2)/(H_v2(a,w_0,w_a))*(np.exp(3*w_a*(a-1.)))/(a**(3.*(1+w_0+w_a)))
def omg_vk(a,w_0,w_a):
    return (1./3.)*(-eta_fxn(a,w_0,w_a)*H0**2)/(H_v2(a,w_0,w_a)*a**2)*(5*omega_m + 2*omega_l*(np.exp(3*w_a*(a-1.)))/(a**(3.*(w_0+w_a))))


#################################################################################
# NEWWWWWWWWWWW CPL PARAMETRIZATION FUNCTIONS
#################################################################################
# BACKGROUND DOESNT CHANGE
# this version we change the a_vs incorrectly assigned in the last bracket to a0 = 0.001
# which means for this version div_time has to be early
# THIS IS THE CORRECT MOST COMPLETE VERSION BUT WE WILL DROP THE OMGL0 TERM TOO
# VOID
a_i = 0.1
def new_da_dt_v(a,t,w_0,w_a):
    eta_value = eta_fxn(a,w_0,w_a)
    deriv = H0*(omega_m/a+omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a))
    -eta_value/3.*(5.*omega_m+2.*omega_l*np.exp(3.*w_a*(a_i-1.))/a_i**(3.*(w_0+w_a))))**(1./2.)
    return deriv
def new_H_v2(a,w_0,w_a):
    val = H0**2*(omega_m/a**3 
    + omega_l*np.exp(3.*w_a*(a-1.))/a**(3*(1+w_0+w_a)) 
    - eta_fxn(a,w_0,w_a)/(3.*a**2)*(5.*omega_m+2.*omega_l*np.exp(3.*w_a*(a_i-1.))/a_i**(3*(w_0+w_a))))
    return val
def new_omg_vM(a,w_0,w_a):
    return (omega_m*H0**2)/(new_H_v2(a,w_0,w_a)*a**3)
def new_omg_vL(a,w_0,w_a):
    return (omega_l*H0**2)/(new_H_v2(a,w_0,w_a))*(np.exp(3*w_a*(a-1.)))/(a**(3.*(1+w_0+w_a)))
def new_omg_vk(a,w_0,w_a):
    return (1./3.)*(-eta_fxn(a,w_0,w_a)*H0**2)/(new_H_v2(a,w_0,w_a)*a**2)*(5*omega_m + 2*omega_l*(np.exp(3*w_a*(a_i-1.)))/(a_i**(3.*(w_0+w_a))))


#################################################################################
# DROP CPL PARAMETRIZATION FUNCTIONS
#################################################################################
# BACKGROUND DOESNT CHANGE
# this version we drop the OmgL0 term at G+V 2004's eq. 9
def drop_da_dt_v(a,t,w_0,w_a):
    eta_value = eta_fxn(a,w_0,w_a)
    deriv = H0*(omega_m/a+omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a))
    -eta_value*5.*omega_m/3.)**(1./2.)
    return deriv
def drop_H_v2(a,w_0,w_a):
    val = H0**2*(omega_m/a**3 
    + omega_l*np.exp(3.*w_a*(a-1.))/a**(3*(1+w_0+w_a)) 
    - eta_fxn(a,w_0,w_a)/(3.*a**2)*(5.*omega_m))
    return val
def drop_omg_vM(a,w_0,w_a):
    return (omega_m*H0**2)/(drop_H_v2(a,w_0,w_a)*a**3)
def drop_omg_vL(a,w_0,w_a):
    return (omega_l*H0**2)/(drop_H_v2(a,w_0,w_a))*(np.exp(3*w_a*(a-1.)))/(a**(3.*(1+w_0+w_a)))
def drop_omg_vk(a,w_0,w_a):
    return (-5.*omega_m*eta_fxn(a,w_0,w_a)*H0**2)/(3*drop_H_v2(a,w_0,w_a)*a**2)




# SHARED PROPERTIES
a0 = np.array([0.001])
t_range_bkgd = np.linspace(0.001,15.,101)
trange_void = np.linspace(div_time, 15.,101)


#################################################################################
# LCDM SOLVING IVPs AND OMEGAS
#################################################################################
# LCDM SOLVE BKGD IVP
lcdm_bkgd_sol = solve_ivp(lcdm_da_dt,[0.001,15.],y0=a0,t_eval=t_range_bkgd) 
lcdm_renorm_bkgd = lcdm_bkgd_sol.y[0]/lcdm_bkgd_sol.y[0][-1]

# LCDM SOLVE VOID IVP
lcdm_av0_idx = np.argmin(abs(lcdm_bkgd_sol.t-div_time))
lcdm_av0 =  lcdm_bkgd_sol.y[0][lcdm_av0_idx]
lcdm_void_sol = solve_ivp(lcdm_da_dt_v,[div_time, 15.],y0=[lcdm_av0], t_eval = trange_void)
lcdm_renorm_void = lcdm_void_sol.y[0]/lcdm_bkgd_sol.y[0][-1]
lcdm_alpha = lcdm_renorm_void[-1]
lcdm_delta_t0 = 1./lcdm_alpha**3-1.

print('LCDM bkgd IVP: '+lcdm_bkgd_sol.message)
print('LCDM void IVP: '+lcdm_void_sol.message)
print('(LCDM) alpha, delta_t0',lcdm_alpha,lcdm_delta_t0)

# Bkgd Omegas
lcdm_omega_Ms = []
lcdm_omega_Ls = []
for value in lcdm_renorm_bkgd:
    lcdm_omega_Ms.append(lcdm_omg_M(value))
    lcdm_omega_Ls.append(lcdm_omg_L(value))

# Void Omegas
lcdm_omega_v_Ms = []
lcdm_omega_v_Ls = []
lcdm_omega_v_ks = []
for value in lcdm_renorm_void:
    lcdm_omega_v_Ms.append(lcdm_omg_vM(value))
    lcdm_omega_v_Ls.append(lcdm_omg_vL(value))
    lcdm_omega_v_ks.append(lcdm_omg_vk(value))

#################################################################################
# CPL PARAM SOLVING IVPs AND OMEGAS
#################################################################################
# CPL SOLVE BKGD IVP
bkgd_sol = solve_ivp(da_dt,[0.001,15.],y0=a0,t_eval=t_range_bkgd, args=(w_0,w_a)) 
renorm_bkgd = bkgd_sol.y[0]/bkgd_sol.y[0][-1]

# CPL SOLVE VOID IVP
av0_idx = np.argmin(abs(bkgd_sol.t-div_time))
av0 =  bkgd_sol.y[0][av0_idx]
void_sol = solve_ivp(da_dt_v,[div_time, 15.],y0=[av0], t_eval = trange_void, args=(w_0, w_a))
renorm_void = void_sol.y[0]/bkgd_sol.y[0][-1]
alpha = renorm_void[-1]
delta_t0 = 1./alpha**3-1.

print('bkgd IVP: '+bkgd_sol.message)
print('void IVP: '+void_sol.message)
print('alpha, delta_t0',alpha,delta_t0)

# Bkgd Omegas
omega_Ms = []
omega_Ls = []
for value in renorm_bkgd:
    omega_Ms.append(omg_M(value,w_0,w_a))
    omega_Ls.append(omg_L(value,w_0,w_a))

# Void Omegas
omega_v_Ms = []
omega_v_Ls = []
omega_v_ks = []
for value in renorm_void:
    omega_v_Ms.append(omg_vM(value,w_0,w_a))
    omega_v_Ls.append(omg_vL(value,w_0,w_a))
    omega_v_ks.append(omg_vk(value,w_0,w_a))

#################################################################################
# NEWWWWWWWWWWWW CPL PARAM SOLVING IVPs AND OMEGAS
#################################################################################
# BKGD IS THE SAME
# NEW CPL SOLVE VOID IVP
new_void_sol = solve_ivp(new_da_dt_v,[div_time,15.],y0=[av0], t_eval=trange_void, args=(w_0,w_a))
new_renorm_void = new_void_sol.y[0]/bkgd_sol.y[0][-1]
new_alpha = new_renorm_void[-1]
new_delta_t0 = 1./new_alpha**3-1.

print('NEW void IVP: '+new_void_sol.message)
print('NEW alpha, NEW delta_t0',new_alpha,new_delta_t0)

# NEW Void Omegas
new_omega_v_Ms = []
new_omega_v_Ls = []
new_omega_v_ks = []
for value in new_renorm_void:
    new_omega_v_Ms.append(new_omg_vM(value,w_0,w_a))
    new_omega_v_Ls.append(new_omg_vL(value,w_0,w_a))
    new_omega_v_ks.append(new_omg_vk(value,w_0,w_a))


#################################################################################
# DROP CPL PARAM SOLVING IVPs AND OMEGAS
#################################################################################
# BKGD IS THE SAME
# DROP CPL SOLVE VOID IVP
drop_void_sol = solve_ivp(drop_da_dt_v,[div_time,15.],y0=[av0], t_eval=trange_void, args=(w_0,w_a))
drop_renorm_void = drop_void_sol.y[0]/bkgd_sol.y[0][-1]
drop_alpha = drop_renorm_void[-1]
drop_delta_t0 = 1./drop_alpha**3-1.

print('DROP void IVP: '+drop_void_sol.message)
print('DROP alpha, DROP delta_t0',drop_alpha,drop_delta_t0)

# NEW Void Omegas
drop_omega_v_Ms = []
drop_omega_v_Ls = []
drop_omega_v_ks = []
for value in drop_renorm_void:
    drop_omega_v_Ms.append(drop_omg_vM(value,w_0,w_a))
    drop_omega_v_Ls.append(drop_omg_vL(value,w_0,w_a))
    drop_omega_v_ks.append(drop_omg_vk(value,w_0,w_a))







#################################################################################
# PLOTTING
#################################################################################
# a(t), a_v(t)
fig, ax = plt.subplots(figsize=(8,6))
ax.set_xlabel(r'Time (Gy)')
ax.set_ylabel(r'a(t)')
ax.plot(lcdm_bkgd_sol.t,lcdm_renorm_bkgd,label='$\Lambda$CDM Background')
ax.plot(lcdm_void_sol.t,lcdm_renorm_void,label='$\Lambda$CDM Void')
ax.plot(bkgd_sol.t,renorm_bkgd,label='Background')
#ax.plot(void_sol.t,renorm_void,label='Void')
#ax.plot(new_void_sol.t,new_renorm_void,label='eq. 3.12 Void')
#ax.plot(drop_void_sol.t,drop_renorm_void,label='eq. 3.13 Void')
ax.plot(drop_void_sol.t,drop_renorm_void,label='Void')
ax.legend()


# Omegas
fig1, ax1 = plt.subplots(figsize=(8,6))
ax1.set_xlabel(r'Time (Gy)')
ax1.set_ylim([-0.1,1.1])
ax1.plot(bkgd_sol.t,omega_Ms,label=r'$\Omega_{M}$')
ax1.plot(bkgd_sol.t,omega_Ls,label=r'$\Omega_{\Lambda}$')
#ax1.plot(void_sol.t,omega_v_Ms,label=r'$\Omega_{v,M}$')
#ax1.plot(void_sol.t,omega_v_Ls,label=r'$\Omega_{v,\Lambda}$')
#ax1.plot(void_sol.t,omega_v_ks,label=r'$\Omega_{v,k}$')

# NEW OMGS
#ax1.plot(new_void_sol.t,new_omega_v_Ms,label=r'NEW $\Omega_{v,M}$')
#ax1.plot(new_void_sol.t,new_omega_v_Ls,label=r'NEW $\Omega_{v,\Lambda}$')
#ax1.plot(new_void_sol.t,new_omega_v_ks,label=r'NEW $\Omega_{v,k}$')

# DROP OMGS
ax1.plot(drop_void_sol.t,drop_omega_v_Ms,label=r'$\Omega_{v,M}$',color='k')
ax1.plot(drop_void_sol.t,drop_omega_v_Ls,label=r'$\Omega_{v,\Lambda}$',color='b')
ax1.plot(drop_void_sol.t,drop_omega_v_ks,label=r'$\Omega_{v,k}$',color='y')

ax1.plot(lcdm_bkgd_sol.t,lcdm_omega_Ms,label=r'$\Lambda$CDM $\Omega_{M}$',linestyle='dashed')
ax1.plot(lcdm_bkgd_sol.t,lcdm_omega_Ls,label=r'$\Lambda$CDM $\Omega_{\Lambda}$',linestyle='dashed')
ax1.plot(lcdm_void_sol.t,lcdm_omega_v_Ms,label=r'$\Lambda$CDM $\Omega_{v,M}$',linestyle='dotted',color='k')
ax1.plot(lcdm_void_sol.t,lcdm_omega_v_Ls,label=r'$\Lambda$CDM $\Omega_{v,\Lambda}$',linestyle='dotted',color='b')
ax1.plot(lcdm_void_sol.t,lcdm_omega_v_ks,label=r'$\Lambda$CDM $\Omega_{v,k}$',linestyle='dotted',color='y')
ax1.legend()


#CHECKING ETA BEHAVIOR
lcdm_eta_list = []
eta_list = []
for value in renorm_bkgd:
    eta_list.append(eta_fxn(value,w_0,w_a))
    lcdm_eta_list.append(lcdm_eta_fxn(value))
fig3, ax3 = plt.subplots(figsize=(7,5))
ax3.set_xlabel(r'Time (Gy)')
ax3.set_ylabel(r'$\eta(t)$')
ax3.plot(bkgd_sol.t,eta_list)
ax3.plot(bkgd_sol.t,lcdm_eta_list,alpha = 0.5)
"""
# DIFF PLOT
fig2, axs = plt.subplots(2,1,sharex=True)
baseline_diff = new_renorm_void/new_renorm_void - 1.
drop_diff = drop_renorm_void/new_renorm_void -1.

per_diff = (new_renorm_void-drop_renorm_void)/(abs(new_renorm_void+drop_renorm_void)* 0.5) * 100.

axs[0].plot(new_void_sol.t,new_renorm_void,label='eq. 3.12 Void',linestyle='dotted')
axs[0].plot(drop_void_sol.t,drop_renorm_void,label='eq. 3.13 Void',alpha = 0.5)
axs[0].plot(bkgd_sol.t,renorm_bkgd,label='Background')
axs[0].set_ylabel(r'a(t)')
axs[0].legend()

#axs[1].plot(new_void_sol.t,baseline_diff,linestyle='dotted')
#axs[1].plot(drop_void_sol.t,drop_diff,color='r')
axs[1].plot(drop_void_sol.t,per_diff,color='k')
axs[1].set_xlabel(r'Time (Gy)')
axs[1].set_ylabel(r'\% Difference')
"""

plt.show()