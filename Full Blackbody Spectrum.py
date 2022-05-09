#!/usr/bin/env python
# coding: utf-8


import numpy as np
import math
import matplotlib.pyplot as plt

# Input Parameters
m1 = 2e9   #mass of primary BH [solar masses]
m2 = 1e9   #mass of secondary BH [solar masses]
e = 0.01   #eccentricity
a = 0.01   #semimajor axis [pc]
nu = 100   #frequency emitted [nm]
i = 90   #inclination [deg]
omega = 45   #argument of periaspe [deg]
n = 2   #number of orbits

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
E_ang = np.linspace(0, 2*math.pi, num=500)   #eccentric anomaly [rad]
E_ang = E_ang.tolist()
m = []   #mean anomaly [rad]
t = []   #time for one period [s]
for j in range (len(E_ang)):
    m.append(E_ang[j] - e * math.sin(E_ang[j]))
    t.append(m[j] / w)
    
# Anomalies and Time for n orbits
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


"""
Entire photon spectrum for a minidisk. Unprimed values refer to the minidisk's reference frame, and primed 
values refer to the observer's reference frame. 
"""
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from decimal import Decimal, getcontext
from mpmath import *
mp.dps=10

# Input Parameters
T = 5000   #disk temperature [K]
d = 1500   #distance to BH [ly]

# Constants
G = 6.67e-11   #gravitational constant [m^3/kg*s^2]
c = 3e8   #speed of light [m/s]
h = 6.626e-34   #Planck's consant [J/Hz]
sigma = 5.67e-8   #Stefan-Boltzmann constant [W/m^2*K^4] 
kb = 1.38e-23   #Boltzmann constant [J/K]

# Convert to SI units
d = d * 9.461e15   #distance to BH [m]

# Wavelength and Frequency
lamda_peak = 2.898e-3 / T   #peak wavelength [m]
lamda = np.linspace(200e-9, 1000e-9, num=500)   #wavelengths in UV spectrum [m]
nu = []   #frequency [s^-1]
for j in range (len(lamda)):
    nu.append(c / lamda[j])
    
# Power Received
rs = a / 3   #radius of disk [m] from models in Minidisks in Binary Black Hole Accretion, Geoffrey Ryan and Andrew McFayden
ar_s = math.pi * rs**2   #area of disk/source [m^2]
Pe = sigma * ar_s * T**4   #power emitted [W]

# Number of Photons and Energy
E_p = []   #energy per photon [J]
for j in range (len(lamda)):
    E_p.append(h * nu[j])
    
z1 = []
z2 = []
z3 = []
z4 = []
B_nu = []   #power per frequency
for j in range (len(lamda)):
    z1.append(2 * h * nu[j]**3)
    z2.append(c**2)
    z3.append(exp(h * nu[j] / (kb * T)))
    z4.append(z3[j] + 1)
    B_nu.append(z1[j] / (z2[j] * z4[j]))
    #B_nu.append(2 * h * nu[j]**3 / (c**2 * (exp(h * nu[j] / (k * T)) - 1)))

E_nu = []   #energy per frequency [Js]
E_bin = []   #energy [J]
N = []   #number of photons
Ndot = []   #photon flux
dt = 1   #time interval [s]
for j in range (len(lamda) - 1):
    nu0 = nu[j]
    nu1 = nu[j+1]
    nuavg = 0.5 * (nu0 + nu1)
    E_nu.append(B_nu[j] * dt)
    E_bin.append(E_nu[j] * (nu0 - nu1))
    N.append(E_bin[j] / h / nuavg)
    Ndot.append(N[j] / (nu0 - nu1) / dt)

# Velocity
vtot = []   #total velocity [m/s]
vpar = []   #photon velocity along LOS [m/s]
for j in range (len(t)):
    vpar.append(v1par[j])   #source moving away from observer must have positive velocity
    vtot.append(-v1[j])
   
# Doppler Shifts
beta = []   #beta ratio
betapar = []  #parallel beta ratio
for j in range (len(t)):
    beta.append(vtot[j] / c)
    betapar.append(vpar[j] / c)
    
gamma = []
for j in range (len(t)):
    gamma.append(1 / math.sqrt(1 - beta[j]**2))

nu_prime = np.zeros((len(t), len(lamda)))   #Doppler shifted frequency [s^-1]
t_prime = np.zeros((len(t), len(lamda)))   #Doppler shifted time [s]
lamda_peakprime = np.zeros((len(t), len(lamda)))   #Doppler shifted peak wavelength [m]
lamda_prime = np.zeros((len(t), len(lamda)))   #Doppler shifted wavelengths [m]
for j in range (len(t)):
    for k in range (len(lamda)):
        nu_prime[j,k] = nu[k] / (gamma[j] * (1 + betapar[j]))
        t_prime[j,k] = t[j] * gamma[j] * (1 + betapar[j])
        lamda_prime[j,k] = c / nu_prime[j,k]

E_pprime = np.zeros((len(t), len(lamda)))   #Doppler shifted energy per photon [J]
for j in range (len(t)):
    for k in range (len(lamda)):
        E_pprime[j,k] = (h * nu_prime[j,k])
  
z1_prime = []
z2_prime = []
z3_prime = []
z4_prime = []
z5_prime = []
B_nuprime = np.zeros((len(t), len(lamda)))   #Doppler shifted radiation rate
for j in range (len(t)):
    for k in range (len(lamda)):
        z1_prime.append(2 * h * nu_prime[j,k]**3)
        z2_prime.append(c**2)
        z3_prime.append(h * nu_prime[j,k] / (kb * T))
        z4_prime.append(exp(z3_prime[k]))
        z5_prime.append(z4_prime[k] + 1)
        B_nuprime[j,k] = (z1_prime[k] / (z2_prime[k] * z5_prime[k]))

E_nuprime = np.zeros((len(t), len(lamda) - 1))
E_binprime = np.zeros((len(t), len(lamda) - 1))
N_prime = np.zeros((len(t), len(lamda) - 1))
Ndot_prime = np.zeros((len(t), len(lamda) - 1))   #photon flux
dt_prime = []  #time interval [s]
for j in range (len(t)):
    dt_prime.append(dt * gamma[j] * (1 + betapar[j]))
    for k in range(len(lamda) - 1):
        nu0 = nu_prime[j,k]
        nu1 = nu_prime[j,k+1]
        nuavg = 0.5 * (nu0 + nu1)
        E_nuprime[j,k] = (B_nuprime[j,k] * dt_prime[j])
        E_binprime[j,k] = (E_nuprime[j,k] * (nu0 - nu1))
        N_prime[j,k] = (E_binprime[j,k] / h / nuavg)
        Ndot_prime[j,k] = (N[k] / (nu0 - nu1) / dt_prime[j])
               
# Flux Received
Fr = Pe / (4 * math.pi * d**2)   #flux received [W/m^2]
    
# Convert Time into Years
t_yr = []   #time [yr]
t_primeyr = np.zeros((len(t), len(lamda)))   #Doppler shifted time [yr]
for j in range (len(t)):
    t_yr.append(t[j] / 3.154e7)
    for k in range (len(lamda)):
        t_primeyr[j,k] = (t_prime[j,k] / 3.154e7)
    
# Convert Wavelength into Nanometers
lamda_peaknm = lamda_peak * 1e9   #peak wavelength [nm]
lamda_peaknmprime = lamda_peakprime * 1e9   #Doppler shifted peak wavelength [nm]
lamda_nm = np.zeros((len(t), len(lamda)))   #wavelengths [nm]
lamda_nmprime = np.zeros((len(t), len(lamda)))   #Doppler shifted wavelengths [nm]
for j in range (len(t)):
    for k in range (len(lamda)):
        lamda_nm[j,k] = (lamda[k] * 1e9)
        lamda_nmprime[j,k] = (lamda_prime[j,k] * 1e9) 
  
