#!/usr/bin/env python
# coding: utf-8
"""
Photon spectrum for a monochromatic minidisk. Unprimed values refer to the minidisk's reference frame, and primed 
values refer to the observer's reference frame.
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# Input Parameters
m1 = 2e8   #mass of primary BH [solar masses]
m2 = 1e8   #mass of secondary BH [solar masses]
e = 0.8   #eccentricity
a = 0.01   #semimajor axis [pc]
nu = 100   #frequency emitted [nm]
i = 90   #inclination [deg]
omega = 45   #argument of periaspe [deg]

# Constants
G = 6.67e-11   #gravitational constant [m^3/kg*s^2]
c = 3e8   #speed of light [m/s]

# Convert to SI Units
m1 = m1 * 1.989e+30   #mass of primary BH [kg]
m2 = m2 * 1.989e+30   #mass of secondary BH [kg]
a = a * 3.086e+16   #semimajor axis [m]
nu = nu * 1e-9   #frequency [m]
i = math.radians(i)   #inclination [rad]
omega = math.radians(omega)   #argument of periapse [rad]

# Mass Ratio and Center of Mass
q = m2 / m1   #mass ratio
mtot = m1 + m2
mew = (m1 * m2) / (m1 + m2)   #reduced mass [kg]
r1_0 = -mew * a / m1   #radial position of primary BH wrt CM at t=0 (negative since @ periapse) [m]
r2_0 = mew * a / m2   #radial position of secondary BH wrt CM at t=0 [m]
cm = ((m1 * r1_0) + (m2 * r2_0)) / mtot   #center of mass (should =0) [m]

# Period and Orbital Frequency
P = math.sqrt((a**3 * 4 * (math.pi)**2) / (G * mtot))   #orbital period of BHs [s]
w = 2 * math.pi / P   #orbital frequency of BHs [s^-1] 

# Anomalies and Time
f_0 = 0   #initial true anomoly [rad] 
E_ang = np.linspace(0, 2*math.pi, num=1000)   #eccentric anomaly [rad]
E_ang = E_ang.tolist()
m = []   #mean anomaly [rad]
t = []   #time for one period [s]
for j in range (len(E_ang)):
    m.append(E_ang[j] - e * math.sin(E_ang[j]))
    t.append(m[j] / w)
    
# Anomalies and Time for n orbits
n = 1   #number of orbits
temp2pi = 2 * math.pi
mtemp = []
Etemp = []
ttemp = []
for k in range (n):
    for j in range (len(m)):
        ttemp.append((m[j] + temp2pi) / w)
    temp2pi += 2 * math.pi
for j in range (n):
    Etemp += E_ang
    mtemp += m
    
E_ang = Etemp
m = mtemp
t = ttemp

# Position and Radial Velocity
p1 = abs(r1_0) * (1 + e * math.cos(f_0))   #semi-latus rectum of primary BH [m]
p2 = abs(r2_0) * (1 + e * math.cos(f_0))   #semi-latus rectum of secondary BH [m]
r1_m = []   #radial distance to primary BH from CM [m]
r2_m = []   #radial distance to secondary BH from CM [m]
for j in range (len(t)):
    r1_m.append((-p1 / (1 + e * math.cos(m[j]))))
    r2_m.append((p2 / (1 + e * math.cos(m[j]))))
    
r1dot = []
r2dot = []
for j in range (len(t)):
    r1dot.append(r1_m[j] * e * w * math.sin(m[j]) / (1 + e * math.cos(m[j])))
    r2dot.append(r2_m[j] * e * w * math.sin(m[j]) / (1 + e * math.cos(m[j])))

x1 = []   #x-position of primary BH [m]
x2 = []   #x-position of secondary BH [m]
y1 = []   #y-position of primary BH [m]
y2 = []   #y-position of secondary BH [m]
v1 = []   #radial velocity of primary BH [m/s]
v2 = []   #radial velocity of secondary BH [m/s]
for j in range (len(t)):
    x1.append(r1_m[j] * math.cos(m[j]))
    x2.append(r2_m[j] * math.cos(m[j]))
    y1.append(r1_m[j] * math.sin(m[j]))
    y2.append(r2_m[j] * math.sin(m[j]))
    v1.append((r1_m[j] * w * math.cos(m[j]) + (r1dot[j] * math.sin(m[j]))) * math.sin(i))
    v2.append((r2_m[j] * w * math.cos(m[j]) + (r2dot[j] * math.sin(m[j]))) * math.sin(i))

# Parallel Velocity (observed velocity offset due to orbital motion projected along LOS)
v1parc = np.zeros((len(t), 3))   #velocity components of primary BH parallel to LOS
v2parc = np.zeros((len(t), 3))   #velocity components of secondary BH parallel to LOS
vparhat = np.array([[math.cos(i) * math.sin(omega)], [math.cos(i) * math.cos(omega)], [math.sin(i)]]) #unit vector in line of sight direction
for j in range (len(t)):
    v1parc[j,0] = v1[j] * vparhat[0]   #x component
    v1parc[j,1] = v1[j] * vparhat[1]   #y component
    v1parc[j,2] = v1[j] * vparhat[2]   #z component
    v2parc[j,0] = v2[j] * vparhat[0]
    v2parc[j,1] = v2[j] * vparhat[1]
    v2parc[j,2] = v2[j] * vparhat[2]

v1par = []   #total velocity of primary BH parallel to LOS
v2par = []   #total velocity of secondary BH parallel to LOS
for j in range (len(t)):
    if (v1parc[j,0]>0) & (v1parc[j,1]>0) & (v1parc[j,2]>0):
        v1par.append(math.sqrt(v1parc[j,0]**2 + v1parc[j,1]**2 + v1parc[j,2]**2))
    else:
        v1par.append(-(math.sqrt(v1parc[j,0]**2 + v1parc[j,1]**2 + v1parc[j,2]**2)))
        
for j in range (len(t)):
    if (v2parc[j,0]>0) & (v2parc[j,1]>0) & (v2parc[j,2]>0):
        v2par.append(math.sqrt(v2parc[j,0]**2 + v2parc[j,1]**2 + v2parc[j,2]**2))
    else:
        v2par.append(-(math.sqrt(v2parc[j,0]**2 + v2parc[j,1]**2 + v2parc[j,2]**2)))      

# Beta and Doppler Factors
beta1 = []   #beta ratio for radial velocity
beta2 = []
beta1par = []   #beta ratio for parallel velocity
beta2par = []
D1 = []   #doppler factor for primary BH
D2 = []   #doppler factor for secondary BH

for j in range (len(t)):
    beta1.append(v1[j] / c)
    beta2.append(v2[j] / c)
    beta1par.append(v1par[j] / c)
    beta2par.append(v2par[j] / c)
    
for j in range (len(t)):
    D1.append(((1 - beta1par[j]) / math.sqrt(1 - beta1[j]**2))**(-1))
    D2.append(((1 - beta2par[j]) / math.sqrt(1 - beta2[j]**2))**(-1))
    
# Observed Frequency and Intensity
alpha = -0.44   #spectral index, used value from paper
K = (m2 / mtot) * w * a * math.sin(i) / math.sqrt(1 - e**2)
nuobs1 = []   #observed frequency of primary BH after Doppler shift [s^-1]
nuobs2 = []   #observed frequency of secondary BH after Doppler shift [s^-1]
Iobs1 = []   #observed intensity of primary BH after Doppler shift [W/m^2]
Iobs2 = []   #observed intensity of secondary BH after Doppler shift [W/m^2]
for j in range (len(t)):
    nuobs1.append(D1[j] * nu)
    nuobs2.append(D2[j] * nu)
    Iobs1.append((nuobs1[j] / nu)**3 * K * (nu)**(alpha))
    Iobs2.append((nuobs2[j] / nu)**3 * K * (nu)**(alpha))
I = K * nu**alpha   #actual intensity [W/m^2]

# Input Parameters
T = 10000   #disk temperature [K]
d = 1500   #distace to BH [ly]
ar_o = 28.27   #area of reciever/observer [m^2] (set as telescope with 3m-radius)

# Constants
G = 6.67e-11   #gravitational constant [m^3/kg*s^2]
c = 3e8   #speed of light [m/s]
h = 6.626e-34   #Planck's consant [J/Hz]
sigma = 5.67e-8   #Stefan-Boltzmann constant [W/m^2*K^4] 
k = 1.38e-23   #Boltzmann constant [J/K]

# Convert to SI Units
d = d * 9.461e15   #distace to BH [m]

# Wavelength and Frequency
lamda_peak = 2.898e-3 / T   #peak wavelength [m]
nu = c / lamda_peak   #frequency [s^-1]   
    
# Number of Photons and Energy
E_p = h * c / lamda_peak   #energy per photon [J]
N = sigma * T**4 / E_p   #number of photons    
E = N * E_p   #total energy of photons [J]
   
Ndot = []   #spectral flux per unit wavelength(dN/dEdt) [W/m]
for j in range (len(t)):
    Ndot.append(N / (E * t[j]))

# Power Received
rs = 3 * a * (1 + np.sqrt(2)) / 2   #radius of disk [m] from models in Minidisks in Binary Black Hole Accretion, Geoffrey Ryan and Andrew McFayden
ar_s = math.pi * rs**2   #area of disk/source [m^2]
Pr = sigma * ar_s * T**4   #power received [W]

# Velocity
vpho = []   #shift in photon velocity along LOS [m/s]
for j in range (len(t)):
    vpho.append(-v1par[j])   #negative since source moving away from observer has positive velocity
    
# Doppler Shifts
beta = []   #beta ratio
for j in range (len(t)):
    beta.append(vpho[j] / c)

T_prime = []   #Doppler shifted temperature [K]
nu_prime = []   #Doppler shifted frequency [s^-1]
t_prime = []   #Doppler shifted time [s]
lamda_prime = []   #Doppler shifted wavelength [m]
for j in range (len(t)):
    T_prime.append(T * ((c - vpho[j]) / (c + vpho[j]))**0.5)
    nu_prime.append(nu * (1 + beta[j])**0.5 / (1 - beta[j])**0.5)
    t_prime.append(t[j] / (1 - beta[j]**2)**0.5)
    lamda_prime.append(lamda_peak * (1 - beta[j]))

E_prime = []   #Doppler shifted radiant energy [J]
for j in range (len(t)):
    E_prime.append(sigma * T_prime[j]**4)
                     
N_prime = []   #Doppler shifted number of photons 
for j in range (len(t)):
    N_prime.append(E_prime[j] / (h * nu_prime[j]))
 
Ndot_prime = []   #Doppler shifted spectral flux per unit wavelength[W/m]
for j in range (len(t)):
    Ndot_prime.append((N_prime[j]) / ((E_prime[j]) * (t_prime[j])))
    
Pr_prime = []   #Doppler shifted power received [W]
for j in range (len(t)):
    Pr_prime.append(sigma * ar_s * T_prime[j]**4)  
    
# Convert Time into Years
t_yr = []   #time [yr]
t_primeyr = []   #Doppler shifted time [yr]
for j in range (len(t)):
    t_yr.append(t[j] / 3.154e7)
    t_primeyr.append(t_prime[j] / 3.154e7)
    
# Convert Wavelength into Nanometers
lamda_nm = lamda_peak * 1e9   #wavelength [nm]
lamda_nmprime = []   #Doppler shifted wavelength [nm]
for j in range (len(t)):
    lamda_nmprime.append(lamda_prime[j] * 1e9) 

"""    
# Array of Unprimed Values For Plotting
lamda_nmar = np.linspace(lamda_nm, lamda_nm, num=len(t))
Pr_ar = np.linspace(Pr, Pr, num=len(t))
    
# Plots 
#BH is apporaching observer at t~8yr
plt.figure(1)
plt.plot(t_yr, Ndot, 'r-', t_primeyr, Ndot_prime, 'k--')
plt.legend(["Source Frame", "Observer Frame"], loc="upper right")
plt.xlabel("Time [yr]")
plt.ylabel("Spectral Flux (Ndot) [W/m]")
plt.title("Photon Flux")

plt.figure(2)
plt.plot(t_yr, Pr_ar, 'r-', t_primeyr, Pr_prime, 'k--')
plt.legend(["Source Frame", "Observer Frame"], loc="upper right")
plt.xlabel("Time [yr]")
plt.ylabel("Power Recieved [W]")
plt.title("Power Recieved by Observer")

plt.figure(3)
plt.plot(t_yr, lamda_nmar, 'r-', t_primeyr, lamda_nmprime, 'k--')
plt.legend(["Source Frame", "Observer Frame"], loc="upper right")
plt.xlabel("Time [yr]")
plt.ylabel("Wavelength [nm]")
plt.title("Doppler Shifted Wavelengths")
"""
