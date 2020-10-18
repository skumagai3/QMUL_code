"""
Goal: use sim.py functionality to see if there is a characteristic void size at a given epoch.

Research question: does including a time-dependent term to the DE EoS impact the
void size distribution?
"""
import os
import numpy as np
import statistics
from random import uniform
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

#cutoff_time = float(input('Enter cutoff time: '))
cutoff_time_range = np.linspace(0.002,14.99,num=200)

# Gen 1k div_times in between 0, 15 Gyr
num_voids = int(input('Number of voids to simulate: '))
#num_voids = 1000
div_time_list = []
i = 0
while i < num_voids:
    div_time_list.append(uniform(0.001,15.))
    i += 1

# Gen 1k initial underdensities in between -0.01 and -0.2 
# (are those good values? does it matter?)
delta_ti_list = []
i = 0
while i < num_voids:
    delta_ti_list.append(uniform(-0.2,-0.01))
    i += 1


# Make eta final values based off the GV2004 definition
a_i = 0.1 # Not sure ab this, fiddle w it
eta_list = []
for item in delta_ti_list:
    eta_list.append(item/a_i)

# What about initial void sizes? What units is that even in? h^-1 Mpc?
# r_v_ti range: try 0.1-5 Mpc h^-1
r_v_ti_list = []
i = 0 
while i < num_voids:
    r_v_ti_list.append(uniform(0.9,6.))
    i += 1

def full_evo(div_time,delta_ti,eta_value,r_v_ti,w_0,w_a,cutoff_time): # Returns R_eff at that point in time
    #Fxns:
    def da_dt(a,t):
        return H0*(omega_m/a + omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a)))**(1./2.)
    def H2(a):
        return H0**2*(omega_m/a**3 + omega_l*np.exp(3.*w_a*(a-1.))/a**(3.*(1.+w_0+w_a)))
    def eta_fxn(a):
        return eta_value * (H0**2)/(H2(a)) * (np.exp(3*w_a*(a-1.)))/(a**(3.*(1.+w_0+w_a)))
    def da_dt_v(a,t):
        deriv = H0*(omega_m/a + omega_l*np.exp(3.*w_a*(a-1.))/a**(1.+3.*(w_0+w_a)) 
        - 5.*eta_fxn(a)*omega_m/3.)
        return deriv

    #BKGD
    a0 = np.array([0.01])
    t_range_bkgd = np.linspace(0.1,15.,101)
    bkgd_sol = solve_ivp(da_dt,[0.001,15.],y0=a0,t_eval=t_range_bkgd) 

    #VOID
    trange_void = np.linspace(div_time, 15.,101)
    av0_idx = np.argmin(abs(bkgd_sol.t-div_time))
    av0 = bkgd_sol.y[0][av0_idx]
    void_sol = solve_ivp(da_dt_v,[div_time, 15.],y0=[av0], t_eval = trange_void)
    renorm_void = void_sol.y[0]/bkgd_sol.y[0][-1]
    #alpha = renorm_void[-1]
    #delta_t0 = 1./alpha**3-1.
    ### R_eff- effective radius. R_eff = r_v(t_0) = alpha*r_v(t_i)/a_i
    #final_R_eff = alpha*r_v_ti/a_i

    cut_off_idx = np.argmin(abs(trange_void-cutoff_time))
    cutoff_a_v = renorm_void[cut_off_idx]

    cut_off_R_eff = cutoff_a_v*r_v_ti/a_i
    
    return cut_off_R_eff


R_eff_list0 = []
R_eff_list1 = []
R_eff_list2 = []
R_eff_list3 = []

i = 0 
while i < num_voids:
    cut_off_radius0 = full_evo(div_time_list[i],delta_ti_list[i],eta_list[i],r_v_ti_list[i],-1.1,-0.2,1.)
    cut_off_radius1 = full_evo(div_time_list[i],delta_ti_list[i],eta_list[i],r_v_ti_list[i],-1.1,-0.2,5.)
    cut_off_radius2 = full_evo(div_time_list[i],delta_ti_list[i],eta_list[i],r_v_ti_list[i],-1.1,-0.2,10.)
    cut_off_radius3 = full_evo(div_time_list[i],delta_ti_list[i],eta_list[i],r_v_ti_list[i],-1.1,-0.2,13.5)

    R_eff_list0.append(cut_off_radius0)
    R_eff_list1.append(cut_off_radius1)
    R_eff_list2.append(cut_off_radius2)
    R_eff_list3.append(cut_off_radius3)
    i += 1

"""
print('Mean: '+str(statistics.mean(R_eff_list0)))
print('Mean: '+str(statistics.mean(R_eff_list1)))
print('Mean: '+str(statistics.mean(R_eff_list2)))
print('Mean: '+str(statistics.mean(R_eff_list3)))

print('Median: '+str(statistics.median(R_eff_list0)))
print('Median: '+str(statistics.median(R_eff_list1)))
print('Median: '+str(statistics.median(R_eff_list2)))
print('Median: '+str(statistics.median(R_eff_list3)))

print('Std dev: '+str(statistics.stdev(R_eff_list0)))
print('Std dev: '+str(statistics.stdev(R_eff_list1)))
print('Std dev: '+str(statistics.stdev(R_eff_list2)))
print('Std dev: '+str(statistics.stdev(R_eff_list3)))
"""

# Cutoff time vs median:
median_list = []
mean_list = []
counter = 0
for value in cutoff_time_range:
    cut_r_eff_list = []
    i = 0
    while i < num_voids:
        cut_rad = full_evo(div_time_list[i],delta_ti_list[i],eta_list[i],r_v_ti_list[i],-1.1,-0.2,value)
        cut_r_eff_list.append(cut_rad)
        i += 1
    median_list.append(statistics.median(cut_r_eff_list))
    mean_list.append(statistics.mean(cut_r_eff_list))
    #print(counter)
    counter += 1







# PLOTTING
fig,ax = plt.subplots(figsize=(8,6))
ax.set_xlabel(r'Cutoff Time (Gy)')
ax.set_ylabel(r'h$^{-1}$ Mpc')
ax.plot(cutoff_time_range,median_list,label=r'Median $R_{eff}$')
ax.plot(cutoff_time_range,mean_list,label=r'Mean $R_{eff}$')
ax.legend(prop={'size': 12})

# Number of bins??? # Rice Rule: k = 2 n^(1/3)
num_bins = int(round(2*num_voids**(1/3)))

fig1,axs = plt.subplots(2,2,figsize=(10,8),sharex=True)
fig1.suptitle('Number of voids simulated: '+str(num_voids))

axs[1,0].set_xlabel(r'$R_{eff}$ ($h^{-1}$ Mpc)')
axs[1,1].set_xlabel(r'$R_{eff}$ ($h^{-1}$ Mpc)')

axs[0,0].hist(R_eff_list0,bins=num_bins,label=r'Cutoff time= 1')
axs[0,1].hist(R_eff_list1,bins=num_bins,label=r'Cutoff time= 5')
axs[1,0].hist(R_eff_list2,bins=num_bins,label=r'Cutoff time= 10')
axs[1,1].hist(R_eff_list3,bins=num_bins,label=r'Cutoff time= 13.5')

axs[0,0].set_title(r'Cutoff time= 1')
axs[0,1].set_title(r'Cutoff time= 5')
axs[1,0].set_title(r'Cutoff time= 10')
axs[1,1].set_title(r'Cutoff time= 13.5')


plt.show()