# -*- coding: utf-8 -*-
"""
Created on Tue May 16 12:57:52 2017

@author: quantummonkey
"""

# Reproducing ADAF Model

import numpy as np
import disk
import matplotlib.pyplot as plt
import scipy.optimize as op

Msun = 1.9891 * 10.**33
yr = 3600. * 24. * 365.25   
G = 6.67259 * 10.**-8                 
cc = 2.99792458 * 10.**10
pc = 3.086 * 10. ** 18  
sigma = 5.669 * 10. ** (-5)
kb = 1.381 * 10. ** (-16)
arad = 7.565 * 10. ** (-15)
mu = 0.6
mH = 1.673 * 10. ** (-24)
kes = 0.4
k0 = 6.4 * 10. ** 22
N = 1e4

gamma = 5./3.
mu = 0.6
fadaf = 1. # ADAF advective parameter
aturb = 0.

def pressure(rho, temp):
    
    return rho* kb * temp/ mu/ mH + 1./3. * arad * temp ** 4.
    
def orbitalspeed(r, Mbh):
    
    return np.sqrt(G * Mbh / r **3. ) * r/(r-disk.Rsch(Mbh))
    
def ffspeed(r, Mbh):
    
    return np.sqrt(G * Mbh/ r)

def adaf_temp(x, rho, pres):

    return pressure(rho, x)- pres
    
def adaf_solvea(x, aturb, visc, fadaf):
    
    return 12.*x*(3.-x) - 24.*(aturb/visc) * (1.-x)*x/(x+5.)/(3.-x)**2. - fadaf    
    
# Define the relative parameters
mbhlist = [1.e8, 1.e9, 1.e10]
macclist = [0.001, 0.01, 0.1] # for proportions
vlist = [0.01,0.03,0.1]

# Figure set up
fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# Main routine

# Below is just for testing the basic thin disk equations

mbh0 = 10.
Mbh = mbh0*Msun
macc0 = float(macclist[1])
Macc = macc0*disk.Medd(Mbh)
visc = float(vlist[2])
dr = disk.Rsch(Mbh) 
r = np.linspace(N*dr, 0, N+1) # grow from outside
    
rtr = N  # outer boundary

Tmp = []
dens = []
h = []
vr = []
pres = []
cs = []
lsp = []

for i in range(len(r)-1):    
        
    vff = ffspeed(r[i], Mbh)
    x = r[i]/dr/rtr
    a = op.fsolve(adaf_solvea, 1., args = (aturb, visc, fadaf))
    csi = np.sqrt(2./5.)* (1.- x**a)**0.5 * vff
    vri = visc * (3.-a)/5. * (1.-x**a) * vff
    omegai = np.sqrt((a+5.)/5.) * x **(a/2.) * vff/ r[i]
    hi = csi/omegai
    rhoi = Macc/(4. * np.pi * r[i] * hi * vri)
    pi = rhoi * csi **2./ gamma
    lspi = omegai * r[i] **2. 
                
                
    vr.append(vri)                      
    h.append(hi)            
    dens.append(rhoi)    
    pres.append(pi)
    cs.append(csi)
    lsp.append(lspi)

    adaf_sol = op.root(adaf_temp, 1e4, args = (dens[i], pres[i],))        
    Tmp.append(adaf_sol.x)
    
    
r = np.delete(r,-1)

ax0.plot(np.array(r)/dr, np.array(cs))
ax1.plot(np.array(r)/dr, np.array(lsp))
ax2.plot(np.array(r)/dr, np.array(dens))
ax3.plot(np.array(r)/dr, np.array(Tmp))
    
    
ax0.set_xlabel(r'$R/R_s$', fontsize = 16)
#ax0.set_ylabel(r'$P \rm [g\ cm^{-1}\ s^{-2}]$', fontsize = 16)
ax0.set_ylabel(r'$c_s [cm/s]$', fontsize = 16)
ax0.set_yscale('log')
ax0.set_xscale('log')

ax1.set_xlabel(r'$R/R_s$', fontsize = 16)
ax1.set_ylabel('l', fontsize = 16)
ax1.set_yscale('log')
ax1.set_xscale('log')

ax2.set_xlabel(r'$R/R_s$', fontsize = 16)
ax2.set_ylabel(r'$\rho [\rm g/cm^3]$', fontsize = 16)
ax2.set_yscale('log')
ax2.set_xscale('log')
    
ax3.set_xlabel(r'$R/R_s$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
ax3.set_xscale('log')

plt.legend(loc='upper right', fontsize=16) 
plt.tight_layout()
plt.savefig('adaf_test.png')
plt.cla()