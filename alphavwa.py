# Goal: create plot of alpha vs w_a, create a colormap
# of alpha vs w_0 and w_a
import numpy as np
import statistics
from random import uniform
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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

# Make w_a continuum
num_points = 20
#w_0_range = np.linspace(-1.2,-0.6,num=num_points) # OG
#w_a_range = np.linspace(-1.6,0.8,num=num_points) # OG
w_0_range = np.linspace(-1.15,-0.6,num=num_points) # Lims from Biswas et al 2010 fig. 3 left panel
w_a_range = np.linspace(-1.,0.8,num=num_points)


def full_evo(div_time,w_0,w_a):
    # Background evolution, need for renormalization
    def da_dt(a,t):
        return H0*(omega_m/a + omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a)))**(1./2.)
    a0 = np.array([0.01])
    t_range_bkgd = np.linspace(0.1,15.,101)
    bkgd_sol = solve_ivp(da_dt,[0.001,15.],y0=a0,t_eval=t_range_bkgd) 
    renorm_bkgd = bkgd_sol.y[0]/bkgd_sol.y[0][-1]


    # Void Evolution
    eta_value = -5.
    trange_void = np.linspace(div_time, 15.,101)
    def H2(a):
        return H0**2*(omega_m/a**3 + omega_l*np.exp(3.*w_a*(a-1.))/a**(3.*(1.+w_0+w_a)))
    def eta_fxn(a):
        return eta_value * (H0**2)/(H2(a)) * (np.exp(3*w_a*(a-1.)))/(a**(3.*(1.+w_0+w_a)))
    def da_dt_v(a,t):
        deriv = H0*(omega_m/a + omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a)) 
        - 5.*eta_fxn(a)*omega_m/3.)
        return deriv
    # CPL SOLVE VOID IVP
    av0_idx = np.argmin(abs(bkgd_sol.t-div_time))
    av0 = bkgd_sol.y[0][av0_idx]
    void_sol = solve_ivp(da_dt_v,[div_time, 15.],y0=[av0], t_eval = trange_void)
    renorm_void = void_sol.y[0]/bkgd_sol.y[0][-1]
    alpha = renorm_void[-1]
    delta_t0 = 1./alpha**3-1.
    return alpha#, delta_t0


div_time_choice = 5.

# COLORMAP
Z_array = []
for w_a in w_a_range:
    row = []
    for w_0 in w_0_range:
        row.append(full_evo(div_time_choice,w_0,w_a))
    Z_array.append(row)

# Try to renormalize?



fig1,ax1=plt.subplots(figsize=(8,6))
ax1.set_xlabel(r'$w_0$')
ax1.set_ylabel(r'$w_a$')

minimum = np.min(Z_array)
maximum = np.max(Z_array)

heatmap = ax1.pcolormesh(w_0_range,w_a_range,Z_array,cmap='jet', norm=colors.LogNorm(vmin=minimum, vmax=maximum))
cb = fig1.colorbar(heatmap,ax=ax1)#,extend='max')
cb.set_label(r'$\alpha$')

plt.show()


"""
w_0_choice = -1.
w_a_choice = 0.



alpha_list = []
delta_t0_list = []

for val in w_a_range:
    # set w_0 = -1 and vary w_a
    alpha_list.append(full_evo(div_time_choice,w_0_choice,w_a=val)[0])
    delta_t0_list.append(full_evo(div_time_choice,w_0_choice,w_a=val)[1])

for val in w_0_range:
    # set w_a = 0 and vary w_a
    alpha_list.append(full_evo(div_time_choice,val,w_a_choice)[0])
    delta_t0_list.append(full_evo(div_time_choice,val,w_a_choice)[1])

# PLOTTING
fig,ax = plt.subplots(figsize=(8,6))
ax.set_xlabel(r'$w_a$')
#ax.set_xlabel(r'$w_0$')
ax.set_ylabel(r'$\alpha$')
ax.set_yscale('log')
ax.plot(w_a_range,alpha_list)
#ax.plot(w_0_range,alpha_list)

fig1,ax1 = plt.subplots(figsize=(8,6))
ax1.set_xlabel(r'$w_a$')
#ax1.set_xlabel(r'$w_0$')
ax1.set_ylabel(r'$\delta(t_o)$')
ax1.plot(w_a_range,delta_t0_list)
#ax1.plot(w_0_range,delta_t0_list)



"""