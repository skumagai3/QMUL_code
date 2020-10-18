"""
Gathering all functions to import into sim.py to analyze void size distribution
version 0 9/2
"""
import numpy as np
from scipy.integrate import solve_ivp

# Constants
G = 6.674e-11           #m^3 kg^-1 s^-2
Mpc = 3.086e19			#km

# bkgd present time params
omega_m = 0.30
omega_l = 0.70
omega_r = 0.00
H0 = 70

# HERE IS WHERE TO SET FINAL_ETA_VALUE FOR sim.py
final_eta_value = -25

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
# CPL PARAMETRIZATION FUNCTIONS
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
