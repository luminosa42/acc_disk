# -*- coding: utf-8 -*-

# Reproducing Thin Disk Model

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
N = 200.

def getpressure(rho, temp):
    
    return rho* kb * temp/ mu/ mH + 1./3. * arad * temp ** 4.
    
def soundspeed(rho, p):
    
    return np.sqrt(p/rho)
    
def orbitalspeed(r, Mbh):
    
    return np.sqrt(G * Mbh / r **3. )
    
def scale_height(cs, vk):
    
    return cs/vk
    
def vel(visc, cs, h, r, dr): # derived from mass & angular momentum conservation

    return visc * cs * h/r / (1- np.sqrt(dr/r))
    
def vel_new(Macc, rho, h, r):
    
    return Macc/ (4. * np.pi * r * rho * h)
    
def kappa(rho, temp):
    
    return kes + k0 * rho * temp **(-3.5)
    
def eden(rho, temp):
    
    return 3./2. * rho* kb * temp/ mu/ mH + arad * temp ** 4

def getm(rho, temp, rout, rin, Mbh):
    
    vk = orbitalspeed(rin, Mbh)
    pres = getpressure(rho, temp)
    return np.pi * (rout **2. - rin **2.)/ vk * np.sqrt(rho * pres)
    
def gete(rho, temp, rout, rin, Mbh):
    
    vk = orbitalspeed(rin, Mbh)
    pres = getpressure(rho, temp)
    cs = soundspeed(rho, pres)
    edens = eden(rho, temp)
    return np.pi * (rout **2. - rin **2.)/ vk * cs * edens
    
def flow(r, dr, visc, Mbh):
    
    A = 2. * np.pi * visc / np.sqrt(G * Mbh)
    return A * r ** (3./2.) /(1- np.sqrt(dr/r))

#Equilibrium setting    
def eloss(Mbh, Macc, r, dr):
    
    return 3./2. * G * Mbh * Macc * dr/ r **2. 

# Individual heating & cooling    
def eheat(visc, cs, rho):
    
    return 3./2. * visc * rho * cs **3. 
    
def ecool(rho, temp, h):
    
    kap = kappa (rho, temp)
    return 32. * sigma * temp **4. / (3.* kap * rho * h )

# This is wrong    
def erad(temp, rout, rin):     
    
    return 2 * np.pi * (rout **2. - rin **2.) * sigma * temp **4.
    
def delta(x):
    
    if x == 0:
        return 1
    else:
        return 0  
      
# Below are the equations to solve      
def eqns(x, consta, constb):
    
    rho, temp = x
    pres = getpressure(rho, temp)
    ui = eden(rho, temp)
    
    return (pres ** (3./2.) * rho ** (-1./2.)-consta, ui * (pres/rho) ** (3./2.)-constb)
    
def eqnset2(x, consta, constb):

    [rho,temp] = x
    pres = getpressure(rho, temp)
    opa = kappa(rho, temp)
    
    return [pres ** (3./2.) * rho ** (-1./2.) -consta, temp **4. /opa/ pres **2.- constb]
    
def eqn_in(x, consta, constb):
    
    [rho,temp] = x
    pres = 1./3. * arad * temp ** 4.
    opa = kes
    return [pres ** (3./2.) * rho ** (-1./2.) -consta, temp **4. /opa/ pres **2.- constb]
    
def eqn_mid(x, consta, constb):
    
    [rho,temp] = x
    return [ temp ** 2. * rho ** (-2.)- consta, temp ** (3./2.) * rho -constb]
    
def eqn_out(x, consta, constb):
    
    [rho, temp] = x
    return[ temp ** 5.5 * rho ** (-3.) - consta, temp ** (3./2.) * rho - constb]    
    
# Define the relative parameters
mbhlist = [1.e8, 1.e9, 1.e10]
macclist = [0.001, 0.01, 0.1] # for proportions
vlist = [0.01,0.03,0.1]
tlist = [1e3, 1e4, 1e5]
rholist = [1e-11, 1e-10, 1e-9]

# Figure set up
fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# Main routine

# Below is just for testing the basic thin disk equations

mbh0 = float(mbhlist[1])
Mbh = mbh0*Msun
macc0 = float(macclist[1])
Macc = macc0*disk.Medd(Mbh)
visc = float(vlist[2])
dr = disk.Risco(Mbh) 
r = np.linspace(N*dr, 0, N+1) # grow from outside
    
tout = tlist[1]
rhoout = rholist[1]

Tmp = []
dens = []
h = []
vr = []
pres = []
Tmp_an = []
dens_an = []
h_an = []
vr_an = []
Tmp_ap = []
dens_ap = []
h_ap = []
vr_ap = []
#    m = []
#    E = []

rab = 150.* (visc * mbh0) ** (2./21.) * (macc0) ** (16./21.) / 3.
rbc = 6300. * (macc0) ** (2./3.) / 3.
print (rab, rbc)

for i in range(len(r)-1):    
        
    vk = orbitalspeed(r[i], Mbh)
    mfac = Macc/(4. * np.pi * visc) * (1- np.sqrt(dr/r[i])) * vk **2. 
    efac = 9. * visc / (8. * arad * cc *vk )
#    print (mfac, efac)
    if (i == 0):    
        sol = op.root(eqnset2, [rhoout, tout], args=(mfac,efac,))
    else:
        sol = op.root(eqnset2, [dens[i-1], Tmp[i-1]], args=(mfac,efac,))

    dens.append(sol.x[0])
    Tmp.append(sol.x[1])
    
    # get the other physical quantities
    pi = getpressure(dens[i], Tmp[i])
    pres.append(pi)                        
    csi = soundspeed(dens[i], pres[i])        
    hi = scale_height(csi, vk)
    h.append(hi)
    #vri = vel(visc, csi, hi, r[i], dr)
    vri = vel_new(Macc, dens[i], hi, r[i])
    vr.append(vri) 
    
     # The analytical forms   

    if (r[i]/dr <= 1 ):
         h_an.append(np.nan)
         dens_an.append(np.nan)
         vr_an.append(np.nan)
         Tmp_an.append(np.nan)
            
         h_ap.append(np.nan)
         dens_ap.append(np.nan)
         vr_ap.append(np.nan)
         Tmp_ap.append(np.nan)
            
    elif (r[i]/dr >1 and r[i]/dr <= rab):
        
        f = 1. - np.sqrt(dr/r[i])
        h_an.append(5.5e4 * mbh0 * macc0 * f)
        dens_an.append(9e-4 * (3.*r[i]/dr)**(3./2.) / (visc * macc0 **2. * f **2. * mbh0))
        vr_an.append(7.6e8 * visc * macc0 **2. * (3.*r[i]/dr)**(-5./2.) *f)
        Tmp_an.append(4.9e7 * (visc * mbh0) ** (-1./4.) * (3.*r[i]/dr)**(-3./8.))
            
        Tmp_ap.append((8. * cc * vk/ (visc * arad * kes)) ** (1./4.))
        const = Macc/(4. * np.pi * visc) * f * vk **2. 
        pi_in = 1./3. * arad * Tmp_ap[i] ** 4. 
        dens_ap.append( pi_in **3.  / const **2. )                         
        csi_in = soundspeed(dens_ap[i], pi_in)        
        hi_in = scale_height(csi_in, vk)
        h_ap.append(hi_in)
        vri_in = vel_new(Macc, dens_ap[i], hi_in, r[i])
        vr_ap.append(vri_in)
        
    elif (r[i]/dr > rab and r[i]/dr <= rbc):
        
        f = 1. - np.sqrt(dr/r[i])
        h_an.append(2.7e3 * visc **(-1./10.) * mbh0 ** (9./10.) * macc0 ** (1./5.) * (3*r[i]/dr)**(21./20.) * f ** (1./5.))
        dens_an.append(8. * visc ** (-7./10.) * mbh0 ** (-7./10.)* macc0 ** (2./5.) *(3*r[i]/dr)**(-33./20.) * f ** (2./5.))
        vr_an.append(1.7e6 * visc ** (4./5.) * mbh0 ** (-1./5.)* macc0 ** (2./5.) *(3*r[i]/dr)**(-2./5.) * f ** (-3./5.))
        Tmp_an.append(2.2e8 * visc ** (-1./5.) * mbh0 ** (-1./5.)* macc0 ** (2./5.) *(3*r[i]/dr)**(-9./10.) * f ** (2./5.))
            

        fac1 = (kb/mu/mH) **2.  * kes * 9. * visc / (8. * arad * cc *vk )    
        fac2 = (kb/mu/mH) **(-3./2.)* Macc/(4. * np.pi * visc) * (1- np.sqrt(dr/r[i])) * vk **2. 
 
        sol_mid = op.root(eqn_mid, [1e-10, 1e4], args=(fac1,fac2,)) 
        dens_ap.append(sol_mid.x[0])
        Tmp_ap.append(sol_mid.x[1])        
        pi_mid = dens_ap[i]* kb * Tmp_ap[i]/ mu/ mH                     
        csi_mid = soundspeed(dens_ap[i], pi_mid)        
        hi_mid = scale_height(csi_mid, vk)
        h_ap.append(hi_mid)
        vri_mid = vel_new(Macc, dens_ap[i], hi_mid, r[i])
        vr_ap.append(vri_mid) 
        
    else:     
                
        f = 1. - np.sqrt(dr/r[i])
        h_an.append(1.5e3 * visc **(-1./10.) * mbh0 ** (9./10.) * macc0 ** (3./20.) * (3*r[i]/dr)**(9./8.) * f ** (3./20.))
        dens_an.append(47. * visc ** (-7./10.) * mbh0 ** (-7./10.)* macc0 ** (11./20.) *(3.*r[i]/dr)**(-15./8.) * f ** (11./20.))
        vr_an.append(5.4e5 * visc **(4./5.) * mbh0 ** (-1./5.) * macc0 ** (3./10.) * (3.*r[i]/dr)**(-1./4.) * f ** (-7./10.))
        Tmp_an.append(6.9e7 * visc ** (-1./5.) * mbh0 ** (-1./5.)* macc0 ** (3./10.) *(3.*r[i]/dr)**(-3./4.) * f ** (3./10.))
        
        fac3 = (kb/mu/mH) **2. * k0 * 9. * visc / (8. * arad * cc *vk ) 
        fac2 = (kb/mu/mH) **(-3./2.)* Macc/(4. * np.pi * visc) * (1.- np.sqrt(dr/r[i])) * vk **2.
        sol_out = op.root(eqn_out, [1e-10, 1e4], args=(fac3,fac2,))
        dens_ap.append(sol_out.x[0])
        Tmp_ap.append(sol_out.x[1])  
        pi_out = dens_ap[i]* kb * Tmp_ap[i]/ mu/ mH                     
        csi_out = soundspeed(dens_ap[i], pi_out)        
        hi_out = scale_height(csi_out, vk)
        h_ap.append(hi_out)
        vri_out = vel_new(Macc, dens_ap[i], hi_out, r[i])
        vr_ap.append(vri_out) 
          
        # solve for rho, T
#        if (i != len(r)-2 ):
            
#            mfac = flow(r[i], dr, visc, Mbh) * pres[i]**(3./2.) * dens[i] ** (-1./2.)+ Macc*delta(i)/ flow(r[i+1], dr, visc, Mbh)
#            efac = flow(r[i], dr, visc, Mbh) * pres[i]**(3./2.) * dens[i] ** (-3./2.) * eden (dens[i], Tmp[i]) + eloss(Mbh, Macc, r[i], dr)/ flow(r[i+1], dr, visc, Mbh)
#            print (mfac, efac)
#            rho, temp = op.fsolve(eqnset2, (dens[i], Tmp[i]), args=(mfac,efac,))
#            print (rho, temp)
#            Tmp.append(temp)
#            dens.append(rho)
            
    
r = np.delete(r,-1)

ax0.plot(np.array(r)/dr, np.array(h)/np.array(r), label = 'Thin Disk')
ax1.plot(np.array(r)/dr, np.array(vr),label = 'Thin Disk')
ax2.plot(np.array(r)/dr, np.array(dens),label = 'Thin Disk')
ax3.plot(np.array(r)/dr, np.array(Tmp),label = 'Thin Disk')

ax0.plot(np.array(r)/dr, np.array(h_an)/np.array(r), linewidth = 2, label = 'Analytic')
ax1.plot(np.array(r)/dr, np.array(vr_an),linewidth = 2,label = 'Analytic')
ax2.plot(np.array(r)/dr, np.array(dens_an),linewidth = 2,label = 'Analytic')
ax3.plot(np.array(r)/dr, np.array(Tmp_an),linewidth = 2,label = 'Analytic')
    
ax0.plot(np.array(r)/dr, np.array(h_ap)/np.array(r), label = 'Approximate')
ax1.plot(np.array(r)/dr, np.array(vr_ap),label = 'Approximate')
ax2.plot(np.array(r)/dr, np.array(dens_ap),label = 'Approximate')
ax3.plot(np.array(r)/dr, np.array(Tmp_ap),label = 'Approximate')
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax0.set_ylabel(r'$P \rm [g\ cm^{-1}\ s^{-2}]$', fontsize = 16)
ax0.set_ylabel('H/R', fontsize = 16)
ax0.set_ylim(1e-5,1e-2)
ax0.set_yscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax1.set_ylabel(r'$v_r [cm/s]$', fontsize = 16)
ax1.set_ylim(1,1e5)
ax1.set_yscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax2.set_ylabel(r'$\rho [\rm g/cm^3]$', fontsize = 16)
ax2.set_ylim(1e-10,1e-4)
ax2.set_yscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')

plt.legend(loc='upper right', fontsize=16) 
plt.tight_layout()
plt.savefig('thin_test.png')
plt.cla()

"""
# Below is with the binning

fig, axes = plt.subplots(figsize=(12, 8), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

for k in range(len(rholist)): 

    Mbh = float(mbhlist[1])*Msun
    Macc = macclist[1]*disk.Medd(Mbh)
    visc = float(vlist[2])
    dr = disk.Risco(Mbh) 
    r = np.linspace(100*dr, 0, 101) # grow from outside
    
    # Set boundary conditions from outside
    tout = tlist[1]
    rhoout = rholist[k]

    Tmp = [tout]
    dens = [rhoout]
    h = []
    vr = []
    pres = []
    m = []
    E = []

    for i in range(len(r)-1):    
            
        # get the other physical quantities
        pi = getpressure(dens[i], Tmp[i])
        pres.append(pi)                        
        csi = soundspeed(dens[i], pres[i])
        vk = orbitalspeed(r[i], Mbh)
        hi = scale_height(csi, vk)
        h.append(hi)
        vri = vel(visc, csi, hi, r[i], dr)
        #vri = vel_new(Macc, dens[i], hi, r[i])
        vr.append(vri) 
        m.append(getm(dens[i],Tmp[i], r[i], r[i+1], Mbh))
        E.append(gete(dens[i],Tmp[i], r[i], r[i+1], Mbh))
            
        # solve for m, E
        if (i != len(r)-2 ):
            
            mfac = flow(r[i], dr, visc, Mbh) * pres[i]**(3./2.) * dens[i] ** (-1./2.)+ Macc*delta(i)/ flow(r[i+1], dr, visc, Mbh)
#            efac = flow(r[i], dr, visc, Mbh) * pres[i]**(3./2.) * dens[i] ** (-3./2.) * eden (dens[i], Tmp[i]) - eheat(visc, csi, dens[i]) + ecool(dens[i], Tmp[i], hi)/ flow(r[i+1], dr, visc, Mbh)
            efac = flow(r[i], dr, visc, Mbh) * pres[i]**(3./2.) * dens[i] ** (-3./2.) * eden (dens[i], Tmp[i]) + eloss(Mbh, Macc, r[i], dr)/ flow(r[i+1], dr, visc, Mbh)
#            print (mfac, efac)
            rho, temp = op.fsolve(eqns, (dens[i], Tmp[i]), args=(mfac,efac,))
            print (rho, temp)            

            Tmp.append(temp)
            dens.append(rho)
            
    
    r = np.delete(r,-1)

    #ax0.plot(np.array(r)/dr, np.array(m)/Msun,label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    #ax1.plot(np.array(r)/dr, np.array(E), label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax0.plot(np.array(r)/dr, np.array(h), label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax1.plot(np.array(r)/dr, np.array(vr),label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax2.plot(np.array(r)/dr, np.array(dens),label = r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    ax3.plot(np.array(r)/dr, np.array(Tmp),label =r'$\rho_{\rm out}= %.1e\ \rm g/cm^3$'%rholist[k])
    
    
ax0.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax0.set_ylabel(r'$\Delta M/M_{\odot}$', fontsize = 16)
ax0.set_ylabel('H [cm]', fontsize = 16)
ax0.set_yscale('log')

ax1.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax1.set_ylabel('E [ergs]', fontsize = 16)
ax1.set_ylabel(r'$v_r [cm/s]$', fontsize = 16)
ax1.set_yscale('log')

ax2.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
#ax2.set_ylabel(r'$n [\rm cm^{-3}]$', fontsize = 16)
ax2.set_ylabel(r'$\rho [\rm g/cm^3]$', fontsize = 16)
ax2.set_yscale('log')
    
ax3.set_xlabel(r'$R/R_{\rm ISCO}$', fontsize = 16)
ax3.set_ylabel('T [K]', fontsize = 16)
ax3.set_yscale('log')

plt.legend(loc='best', fontsize=16) 
plt.tight_layout()
plt.savefig('thin2.png')
plt.cla()

"""