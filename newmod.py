# -*- coding: utf-8 -*-

# Functions for radius-dependent un-steady accretions

import numpy as np
import disk
import matplotlib.pyplot as plt
#import sympy as sp

Msun = 1.9891 * 10.**33
yr = 3600. * 24. * 365.25   
G = 6.67259 * 10.**-8                 
cc = 2.99792458 * 10.**10
pc = 3.086 * 10. ** 18  
sigma = 5.669 * 10. ** (-5)
kb = 1.381 * 10. ** (-16)
mu = 0.618
mH = 1.673 * 10. ** (-24)
N = 1000 # How many zones

#-------------------------------------------------------------------------------
#---------------------- GENERAL FUNCTIONS FOR ALL DISKS ------------------------
#-------------------------------------------------------------------------------
def soundspeed(e,m):
    
    return np.sqrt(10./9. * e/m)

def surf_dens(m, rout, rin):
    
    return m/ np.pi/ (rout **2. - rin **2.)
    
def vel(visc, cs, h, r, dr):

    return visc * cs * h/r / (1- np.sqrt(dr/r))

def scale_height(Mbh, cs, r) :

    return np.sqrt(3./5.)* cs * np.sqrt(r**3./G / Mbh)   
    
def flux(r, vr, dens):
    
    return 2 * np.pi * r * vr * dens
    
def temp(e,m):
    
    return 2./3. * e/m * mu * mH/kb
    
def erad(T, rout, rin):     
    
    return 2 * np.pi * (rout **2. - rin **2.) * sigma * T **4.
    
def eloss(Mbh, Macc, r, dr):
    
    return 3./2. * G * Mbh * Macc * dr/ r **2.
    
# some specific terms
def terma(r, dr, visc, Mbh):
    
    A = 20./9. * visc /np.sqrt(G * Mbh)
    return A * r ** (3./2.) /(2*r - dr)/ dr /(1- np.sqrt(dr/r)) 
    
def edens(rho, temp):
    
    return 3./2. * rho * kb * temp / mu / mH
    
def mass_from_temp(temp, e):
    
    return 2./3. * e/temp * mu * mH/kb
    
def cs_from_temp(temp):
    
    return np.sqrt(5./3. * kb * temp/ mu/ mH)
    
def delta(x):
    
    if x == 0:
        return 1
    else:
        return 0
    
#-------------------------------------------------------------------------------
#-------------------------------- MAIN PROGRAM ---------------------------------
#-------------------------------------------------------------------------------

# Define the relative parameters
mbhlist = [1.e8, 1.e9, 1.e10]
macclist = [0.001, 0.01, 0.1] # for proportions
vlist = [0.01,0.03,0.1]
tlist = [1e3, 1e4, 1e5]
#elist = [0.01, 0.1, 1, 10, 100]
#rholist = [1e-25, 1e-20, 1e-15]
rholist = [1e-15, 1e-10, 1e-5]
"""
# Figure set up
fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# Main routine
for k in range(len(mbhlist)): 

    Mbh = float(mbhlist[k])*Msun
#    Macc = float(macclist[1])*Msun/yr
    Macc = macclist[1]*disk.Medd(Mbh)
#    print(Macc*yr/Msun)
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(100*dr, 0, 101) # grow from outside
    
    # set the initial m0, E0
    #mdisk = disk.ssd_diskmass(visc, Macc, Mbh, 5./27., 2./3.)
    #m = [mdisk]
    #E = [1.e60]

    # Set boundary conditions from outside
    tout = tlist[1]
    rhoout = rholist[1]
    csout = cs_from_temp(tout)
    hout = scale_height(Mbh, csout, 100*dr)
    eout = edens(rhoout, tout) * np.pi * hout * dr **2. *100  
    mout = mass_from_temp(tout, eout)
    print (eout, mout/Msun)

    Tmp = [tout]
    m = [mout]
    E = [eout]
    h = []
    vr = []
    dens = []
    eden = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities                        
        rhoi = surf_dens(m[i], r[i], r[i+1])
        dens.append(rhoi)    
        csi = soundspeed(E[i], m[i])
        hi = scale_height(Mbh, csi, r[i])
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        vr.append(vri)    
        ui = surf_dens(E[i], r[i], r[i+1])
        eden.append(ui)
    
        # solve for E_i+1 and m_i+1
        if (i != len(r)-2 ):
            newE = (terma(r[i], dr, visc, Mbh) * E[i] + Macc*delta(i))/terma(r[i+1], dr, visc, Mbh)
            E.append(newE)
#        newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + erad(Tmp[i], r[i], r[i+1]))
            newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + eloss(Mbh, Macc, r[i], dr))
            print (eloss(Mbh, Macc, r[i], dr), newE)
            m.append(newm)
            Tmp.append(temp(newE,newm))
    
    r = np.delete(r,-1)

    ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$M_{\rm BH} = %.1e M_{\odot}$'%mbhlist[k])
    ax1.plot(np.array(r)/dr, np.array(E), label = r'$M_{\rm BH} = %.1e M_{\odot}$'%mbhlist[k])
    ax2.plot(np.array(r)/dr, np.array(dens)/np.array(h)/mu/mH, label = r'$M_{\rm BH} = %.1e M_{\odot}$'%mbhlist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp), label = r'$M_{\rm BH} = %.1e M_{\odot}$'%mbhlist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_yscale('log')
#ax0.set_xscale('log')
#ax0.set_xlim(1e-2,1e2)

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_yscale('log')
#ax1.set_xscale('log')
#ax1.set_xlim(1e-2,1e2)

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax2.set_ylabel(r'$\rho [\rm g/cm^3]$', fontsize = 16)
ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_yscale('log')
#ax2.set_xscale('log')
#ax2.set_xlim(1e-2,1e2)
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
#ax3.set_xscale('log')
#ax3.set_xlim(1e-2,1e2)

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('mbh.png')
plt.cla()
# Evolution term
#for i in range(len(r)-1):    
#    mterm = flux(r[i+1], vr[i+1], dens[i+1]) - flux(r[i], vr[i], dens[i])
#    eterm = flux(r[i+1], vr[i+1], eden[i+1]) - flux(r[i], vr[i], eden[i]) - Erad(Tmp[i], r[i+1], r[i])

fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

for k in range(len(macclist)): 

    Mbh = float(mbhlist[1])*Msun
    Macc = macclist[k]*disk.Medd(Mbh)
#    Macc = float(macclist[k])*Msun/yr
#    print(Macc*yr/Msun)
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(100*dr, 0, 101) # grow from outside

    # Set boundary conditions from outside
    tout = tlist[1]
    rhoout = rholist[1]
    csout = cs_from_temp(tout)
    hout = scale_height(Mbh, csout, 100*dr)
    eout = edens(rhoout, tout) * np.pi * hout * dr **2. *100  
    mout = mass_from_temp(tout, eout)
    
    Tmp = [tout]
    m = [mout]
    E = [eout]
    h = []
    vr = []
    dens = []
    eden = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities                        
        rhoi = surf_dens(m[i], r[i], r[i+1])
        dens.append(rhoi)    
        csi = soundspeed(E[i], m[i])
        hi = scale_height(Mbh, csi, r[i])
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        vr.append(vri)    
        ui = surf_dens(E[i], r[i], r[i+1])
        eden.append(ui)
    
        # solve for E_i+1 and m_i+1
        if (i != len(r)-2 ):
            newE = (terma(r[i], dr, visc, Mbh) * E[i] + Macc*delta(i))/terma(r[i+1], dr, visc, Mbh)
            E.append(newE)
#        newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + erad(Tmp[i], r[i], r[i+1]))
            newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + eloss(Mbh, Macc, r[i], dr))
            m.append(newm)
            Tmp.append(temp(newE,newm))
    
    r = np.delete(r,-1)

    ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$M_{\rm BH} = %.3f \dot{M}_{\rm Edd}$'%macclist[k])
    ax1.plot(np.array(r)/dr, np.array(E), label = r'$M_{\rm BH} = %.3f \dot{M}_{\rm Edd}$'%macclist[k])
    ax2.plot(np.array(r)/dr, np.array(dens)/np.array(h)/mu/mH, label = r'$M_{\rm BH} = %.3f \dot{M}_{\rm Edd}$'%macclist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp), label = r'$M_{\rm BH} = %.3f \dot{M}_{\rm Edd}$'%macclist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_yscale('log')
#ax0.set_xscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_yscale('log')
#ax1.set_xscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax2.set_ylabel(r'$\rho [\rm g/cm^3]$', fontsize = 16)
ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_yscale('log')
#ax2.set_xscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
#ax3.set_xscale('log')

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('macc.png')
plt.cla()   

fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

for k in range(len(vlist)): 

    Mbh = float(mbhlist[1])*Msun
    Macc = macclist[1]*disk.Medd(Mbh)  
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(100*dr, 0, 101) # grow from outside

    # Set boundary conditions from outside
    tout = tlist[1]
    rhoout = rholist[1]
    csout = cs_from_temp(tout)
    hout = scale_height(Mbh, csout, 100*dr)
    eout = edens(rhoout, tout) * np.pi * hout * dr **2. *100  
    mout = mass_from_temp(tout, eout)
    #print (eout, mout/Msun)

    Tmp = [tout]
    m = [mout]
    E = [eout]
    h = []
    vr = []
    dens = []
    eden = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities                        
        rhoi = surf_dens(m[i], r[i], r[i+1])
        dens.append(rhoi)    
        csi = soundspeed(E[i], m[i])
        hi = scale_height(Mbh, csi, r[i])
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        vr.append(vri)    
        ui = surf_dens(E[i], r[i], r[i+1])
        eden.append(ui)
    
        # solve for E_i+1 and m_i+1
        if (i != len(r)-2 ):
            newE = (terma(r[i], dr, visc, Mbh) * E[i] + Macc*delta(i))/terma(r[i+1], dr, visc, Mbh)
            E.append(newE)
#        newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + erad(Tmp[i], r[i], r[i+1]))
            newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + eloss(Mbh, Macc, r[i], dr))
            m.append(newm)
            Tmp.append(temp(newE,newm))
    
    r = np.delete(r,-1)

    ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$\alpha = %.2f$'%vlist[k])
    ax1.plot(np.array(r)/dr, np.array(E), label = r'$\alpha = %.2f$'%vlist[k])
    ax2.plot(np.array(r)/dr, np.array(dens)/np.array(h)/mu/mH, label = r'$\alpha = %.2f$'%vlist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp), label = r'$\alpha = %.2f$'%vlist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_yscale('log')
#ax0.set_xscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_yscale('log')
#ax1.set_xscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_yscale('log')
#ax2.set_xscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
#ax3.set_xscale('log')

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('visc.png')
plt.cla() 
"""
fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

for k in range(len(tlist)): 

    Mbh = float(mbhlist[1])*Msun
    Macc = macclist[1]*disk.Medd(Mbh)  
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(N*dr, 0, N+1) # grow from outside

    # Set boundary conditions from outside
    tout = tlist[k]
    rhoout = rholist[1]
#    efac = G * Mbh * Msun/100/dr
    csout = cs_from_temp(tout)
    hout = scale_height(Mbh, csout, N*dr) 
    eout = edens(rhoout, tout) * np.pi * hout * dr **2. *N 
    mout = mass_from_temp(tout, eout)
    print (eout, mout/Msun)

    Tmp = [tout]
    m = [mout]
    E = [eout]
    h = []
    vr = []
    dens = []
    eden = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities                        
        rhoi = surf_dens(m[i], r[i], r[i+1])
        dens.append(rhoi)    
        csi = soundspeed(E[i], m[i])
        hi = scale_height(Mbh, csi, r[i])
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        vr.append(vri)    
        ui = surf_dens(E[i], r[i], r[i+1])
        eden.append(ui)
    
        # solve for E_i+1 and m_i+1
        if (i != len(r)-2 ):
            newE = (terma(r[i], dr, visc, Mbh) * E[i] + Macc*delta(i))/terma(r[i+1], dr, visc, Mbh)
            E.append(newE)
#        newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + erad(Tmp[i], r[i], r[i+1]))
            newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + eloss(Mbh, Macc, r[i], dr))
            m.append(newm)
            Tmp.append(temp(newE,newm))
    
    r = np.delete(r,-1)

    ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$T_{\rm out}= %.1e\ K$'%tlist[k])
    ax1.plot(np.array(r)/dr, np.array(E), label = r'$T_{\rm out}= %.1e\ K$'%tlist[k])
    ax2.plot(np.array(r)/dr, np.array(dens)/np.array(h)/mu/mH, label = r'$T_{\rm out}= %.1e\ K$'%tlist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp), label = r'$T_{\rm out}= %.1e\ K$'%tlist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_yscale('log')
#ax0.set_xscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_yscale('log')
#ax1.set_xscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_yscale('log')
#ax2.set_xscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
#ax3.set_xscale('log')

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('bound1.png')
plt.cla()     

fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

for k in range(len(rholist)): 

    Mbh = float(mbhlist[1])*Msun
    Macc = macclist[1]*disk.Medd(Mbh)  
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(N*dr, 0, N+1) # grow from outside

    # Set boundary conditions from outside
    tout = tlist[1]
    rhoout = rholist[k]
    csout = cs_from_temp(tout)
    hout = scale_height(Mbh, csout, N*dr)
    eout = edens(rhoout, tout) * np.pi * hout * dr **2. *N  
    mout = mass_from_temp(tout, eout)
    print (eout, mout/Msun)

    Tmp = [tout]
    m = [mout]
    E = [eout]
    h = []
    vr = []
    dens = []
    eden = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities                        
        rhoi = surf_dens(m[i], r[i], r[i+1])
        dens.append(rhoi)    
        csi = soundspeed(E[i], m[i])
        hi = scale_height(Mbh, csi, r[i])
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        vr.append(vri)    
        ui = surf_dens(E[i], r[i], r[i+1])
        eden.append(ui)
    
        # solve for E_i+1 and m_i+1
        if (i != len(r)-2 ):
            newE = (terma(r[i], dr, visc, Mbh) * E[i] + Macc*delta(i))/terma(r[i+1], dr, visc, Mbh)
            E.append(newE)
#        newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + erad(Tmp[i], r[i], r[i+1]))
            newm = newE **2. * terma(r[i+1], dr, visc, Mbh)/( terma(r[i], dr, visc, Mbh) * E[i] **2. /m[i] + eloss(Mbh, Macc, r[i], dr))
            m.append(newm)
            Tmp.append(temp(newE,newm))
    
    r = np.delete(r,-1)

    ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax1.plot(np.array(r)/dr, np.array(E), label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax2.plot(np.array(r)/dr, np.array(dens)/np.array(h)/mu/mH,label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp),label =r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_yscale('log')
#ax0.set_xscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_yscale('log')
#ax1.set_xscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_yscale('log')
#ax2.set_xscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')
#ax3.set_xscale('log')

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('bound2.png')
plt.cla()     