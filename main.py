# -*- coding: utf-8 -*-


# Auxilary sink particle program
import numpy as np
import matplotlib.pyplot as plt
import disk    # functions from different disk models

#-------------------------------------------------------------------------------
#-------------------------- PHYSICAL CONSTANTS ---------------------------------
#-------------------------------------------------------------------------------

G = 6.67259 * 10.**-8                 # Units cm^3 g^-1 s^-2
cc = 2.99792458 * 10.**10              # Units cm s^-1
Msun = 1.9891 * 10.**33               # Units g
yr = 3600. * 24. * 365.25             # Units s

# Estimate for viscosities - unitless
visc_SSDob = 0.1                      # Viscosity for thin disk outburst state
visc_SSDqu = 0.01                     # Viscosity for thin disk quiescent state
visc_ADAF  = 0.25                     # Viscosity for an ADAF

#-------------------------------------------------------------------------------
#-------------------------- GENERAL FUNCTIONS ---------------------------------
#-------------------------------------------------------------------------------

#   Get accretion rate for a specific time 
def getaccrate_bondi(Mbh, rho, cs):

#    Bondi-alpha model    
    alpha = 100
    macc = alpha * 4.*np.pi*G**2.* Mbh**2. * rho / cs**3.    
    return macc
    
def getaccrate_period(macc0, P, amp, time):
    
    return macc0 + amp* macc0* np.sin(2* np.pi * time / P )

def getaccrate_linear(macc0, time, tfac):

    return macc0*(1+time/tfac) ** 10.
    
def getaccrate_exp(macc0, time, tfac):
    
    return macc0*np.exp(time/tfac)    
    
#   Constraint on timestep
def timestep(dt, visc, Mbh):

        
    return dt
    
#   Compute disk accretion rate
#   This part is dependent on specific disk model     
def diskacc(visc, mod, macc_new, macc_old, Mbh, dM_jet, ind1, ind2, dt):
    
    if (mod == 'ssd'):
        dM_disk = disk.ssd_diskacc(visc, macc_new, macc_old, dt, Mbh, dM_jet, ind1, ind2 )
    elif (mod == 'adaf'):
        dM_disk = disk.adaf_diskacc(visc, dM_acc, Mbh)
    else:
        print("Invalid disk model, no accretion is calculated")
        dM_disk = 0.
    
    return dM_disk
    

#   Disk evolution
def diskevol(macc_new, macc_old, dM_jet, Mbh, Mdisk, dt, visc, mod, ind1, ind2):
    
    dM_disk = diskacc(visc, mod, macc_new, macc_old, Mbh, dM_jet, ind1, ind2, dt)
    dM_bh = macc_new - dM_disk - dM_jet
    if dM_bh < 0 :
        print('Black hole losing mass!', dM_bh/Msun*yr)
        
    Mbh = Mbh + dM_bh*dt
    Mdisk = Mdisk + dM_disk*dt

    return dM_disk, dM_bh, Mdisk, Mbh

#   Main Program

# These are the parameters to vary
vlist = [visc_SSDqu, visc_SSDob, visc_ADAF]
mbhlist = [1.e8, 1.e9, 1.e10]
macclist = [0.1, 1., 10.]
ex1 = [1./27., 1./9.,5./27.]
ex2 = [1./3., 1./2.,2./3.]

#   Custom Initialization
mod = 'ssd'
dM_jet = 0.
dt = 1e5*yr   # initial dt
rho = 1.e-26
cs = 2.e8
P = 1e7*yr
amp = 0.1
#tloss = 1e9*yr

# Figure set up
fig, axes = plt.subplots(figsize=(16, 12), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# choose parameter to vary here
for k in range(len(vlist)): 
    visc = vlist[k]
    mbh_init = float(mbhlist[1])*Msun
    #macc_init = getaccrate(mbh_init, rho, cs)
    macc_init = float(macclist[1])*Msun/yr
    ind1 = ex1[2]
    ind2 = ex2[2]
    #   General Initialization
    times = [0.0] 
    steps = 1000   
    dM_bh = [0.0]
    dM_disk = [0.0]    
    dM_acc = [macc_init]
    mdisk_init = disk.ssd_diskmass(visc, macc_init, mbh_init, ind1, ind2)
    Mbh = [mbh_init]
    Mdisk = [mdisk_init]

#   Evolution Started!
    for i in range(steps):    
    
        macc = getaccrate_exp(macc_init, times[i], P)
#        macc = getaccrate_period(macc_init,P, amp, times[i])
#        macc = macc_init
        dM_acc.append(macc)  
        a, b, c, d = diskevol(dM_acc[i+1], dM_acc[i], dM_jet, Mbh[i], Mdisk[i], dt, visc, mod, ind1, ind2)
#        b, a, d, c = disk.SSD_MEvol( Mbh[i], Mdisk[i], dM_acc[i+1], dM_acc[i], dM_jet, dt, visc, ind1, ind2)
    
        times.append( times[-1] + dt )
        dM_disk.append(a)
        dM_bh.append(b) 
        Mdisk.append(c)
        Mbh.append(d)
    
#    mod = testregion(visc,b,d)    
    # Get dt for the next step (optional)
    #dt = timestep()
       
#   Output and Plotting (for one specific case) 
       
    ax0.plot(np.array(times)/yr, (np.array(Mbh)-mbh_init)/Msun, label = r'$\alpha = %.2f$'%vlist[k])
    ax1.plot(np.array(times)/yr, (np.array(Mdisk)-mdisk_init)/Msun, label = r'$\alpha = %.2f$'%vlist[k])
    ax2.plot(np.array(times)/yr, np.array(dM_bh)*yr/Msun, label = r'$\alpha = %.2f$'%vlist[k])    
    ax3.plot(np.array(times)/yr, np.array(dM_disk)*yr/Msun, label = r'$\alpha = %.2f$'%vlist[k])
 
 
# Modification of plots etc.
ax0.set_yscale('log')
ax0.set_xlabel('t [yr]', fontsize=16)
ax0.set_ylabel(r'$\Delta M_{BH} [M_{\odot}]$', fontsize = 16)
ax0.set_title('Black Hole Mass - Initial Mass',fontsize=16)

ax1.set_yscale('log')
ax1.set_xlabel('t [yr]' ,fontsize=16)
ax1.set_ylabel(r'$\Delta M_{\rm disk} [M_{\odot}]$', fontsize=16)
ax1.set_title('Disk Mass - Initial Mass',fontsize=16)    

#ax2.set_ylim(1-0.0013,1-0.0008)
ax2.set_xlabel('t [yr]',fontsize=16)
ax2.set_ylabel(r'$\dot{M}_{BH} [M_{\odot}/yr]$',fontsize=16)
ax2.set_title('Black Hole Accretion Rate',fontsize=16)

ax3.set_xlabel('t [yr]',fontsize=16)
ax3.set_ylabel(r'$\dot{M}_{\rm disk} [M_{\odot}/yr]$',fontsize=16)
#ax3.set_ylim(0.0008,0.0013)
ax3.set_title('Disk Accretion Rate',fontsize=16)

plt.legend(loc='best', fontsize=16)    
plt.tight_layout()
plt.savefig('prop_all_visc.png')
plt.cla()

# Figure set up
fig, axes = plt.subplots(figsize=(16, 12), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# choose parameter to vary here
for k in range(len(mbhlist)): 
    visc = vlist[1]
    mbh_init = float(mbhlist[k])*Msun
    macc_init =float(macclist[1])*Msun/yr
    ind1 = ex1[2]
    ind2 = ex2[2]
    #   General Initialization
    times = [0.0] 
    steps = 1000   
    dM_bh = [0.0]
    dM_disk = [0.0]    
    dM_acc = [macc_init]
    mdisk_init = disk.ssd_diskmass(visc, macc_init, mbh_init, ind1, ind2)
    Mbh = [mbh_init]
    Mdisk = [mdisk_init]

    for i in range(steps):    
    
#        macc = macc_init
#        macc = getaccrate_period(macc_init,P, amp, times[i])
        macc = getaccrate_exp(macc_init, times[i], P)
        dM_acc.append(macc)  
        a, b, c, d = diskevol(dM_acc[i+1], dM_acc[i], dM_jet, Mbh[i], Mdisk[i], dt, visc, mod, ind1, ind2)
#        b, a, d, c = disk.SSD_MEvol( Mbh[i], Mdisk[i],dM_acc[i+1], dM_acc[i], dM_jet, dt, visc, ind1, ind2)
    
        times.append( times[-1] + dt )
        dM_disk.append(a)
        dM_bh.append(b) 
        Mdisk.append(c)
        Mbh.append(d)
    
#    mod = testregion(visc,b,d)    
    # Get dt for the next step (optional)
    #dt = timestep()
       
#   Output and Plotting (for one specific case) 
       
    ax0.plot(np.array(times)/yr, (np.array(Mbh)-mbh_init)/Msun, label = r'$M_{\rm BH} = %e M_{\odot}$'%mbhlist[k])
    ax1.plot(np.array(times)/yr, (np.array(Mdisk)-mdisk_init)/Msun, label = r'$M_{\rm BH} = %e M_{\odot}$'%mbhlist[k])
    ax2.plot(np.array(times)/yr, np.array(dM_bh)*yr/Msun, label = r'$M_{\rm BH} = %e M_{\odot}$'%mbhlist[k])    
    ax3.plot(np.array(times)/yr, np.array(dM_disk)*yr/Msun, label = r'$M_{\rm BH} = %e M_{\odot}$'%mbhlist[k])
 
 
# Modification of plots etc.
ax0.set_yscale('log')
ax0.set_xlabel('t [yr]', fontsize=16)
ax0.set_ylabel(r'$\Delta M_{BH} [M_{\odot}]$', fontsize = 16)
ax0.set_title('Black Hole Mass - Initial Mass',fontsize=16)

ax1.set_yscale('log')
ax1.set_xlabel('t [yr]' ,fontsize=16)
ax1.set_ylabel(r'$\Delta M_{\rm disk} [M_{\odot}]$', fontsize=16)
ax1.set_title('Disk Mass - Initial Mass',fontsize=16)    

#ax2.set_ylim(1-0.0025,1)
ax2.set_xlabel('t [yr]',fontsize=16)
ax2.set_ylabel(r'$\dot{M}_{BH} [M_{\odot}/yr]$',fontsize=16)
ax2.set_title('Black Hole Accretion Rate',fontsize=16)

ax3.set_xlabel('t [yr]',fontsize=16)
ax3.set_ylabel(r'$\dot{M}_{\rm disk} [M_{\odot}/yr]$',fontsize=16)
ax3.set_title('Disk Accretion Rate',fontsize=16)

plt.legend(loc='best', fontsize=16)    
plt.tight_layout()
plt.savefig('prop_all_mbh.png')
plt.cla()

# Figure set up
fig, axes = plt.subplots(figsize=(16, 12), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# choose parameter to vary here
for k in range(len(macclist)): 
    visc = vlist[1]
    mbh_init = float(mbhlist[1])*Msun
    macc_init = float(macclist[k])*Msun/yr
    ind1 = ex1[2]
    ind2 = ex2[2]
    #   General Initialization
    times = [0.0] 
    steps = 1000   
    dM_bh = [0.0]
    dM_disk = [0.0]    
    dM_acc = [macc_init]
    mdisk_init = disk.ssd_diskmass(visc, macc_init, mbh_init, ind1, ind2)
    Mbh = [mbh_init]
    Mdisk = [mdisk_init]

#   Evolution Started!
    for i in range(steps):    
    
#        macc = macc_init
#        macc = getaccrate_period(macc_init,P, amp, times[i])
#        macc = getaccrate_linear(macc_init, times[i], tloss)
        macc = getaccrate_exp(macc_init, times[i], P)
        dM_acc.append(macc)  
        a, b, c, d = diskevol(dM_acc[i+1], dM_acc[i], dM_jet, Mbh[i], Mdisk[i], dt, visc, mod, ind1, ind2)
 #       b, a, d, c = disk.SSD_MEvol( Mbh[i], Mdisk[i], dM_acc[i+1], dM_acc[i], dM_jet, dt, visc, ind1, ind2)
    
        times.append( times[-1] + dt )
        dM_disk.append(a)
        dM_bh.append(b) 
        Mdisk.append(c)
        Mbh.append(d)
    
#    mod = testregion(visc,b,d)    
    # Get dt for the next step (optional)
    #dt = timestep()
       
#   Output and Plotting (for one specific case) 
       
    ax0.plot(np.array(times)/yr, (np.array(Mbh)-mbh_init)/Msun, label = r'$\dot{M}_{\rm acc} = %.1f M_{\odot}/yr$'%macclist[k])
    ax1.plot(np.array(times)/yr, (np.array(Mdisk)-mdisk_init)/Msun, label = r'$\dot{M}_{\rm acc} = %.1f M_{\odot}/yr$'%macclist[k])
    ax2.plot(np.array(times)/yr, np.array(dM_bh)*yr/Msun, label = r'$\dot{M}_{\rm acc} = %.1f M_{\odot}/yr$'%macclist[k])    
    ax3.plot(np.array(times)/yr, np.array(dM_disk)*yr/Msun, label = r'$\dot{M}_{\rm acc} = %.1f M_{\odot}/yr$'%macclist[k])
 
 
# Modification of plots etc.
ax0.set_yscale('log')
ax0.set_xlabel('t [yr]', fontsize=16)
ax0.set_ylabel(r'$\Delta M_{BH} [M_{\odot}]$', fontsize = 16)
ax0.set_title('Black Hole Mass - Initial Mass',fontsize=16)

ax1.set_yscale('log')
ax1.set_xlabel('t [yr]' ,fontsize=16)
ax1.set_ylabel(r'$\Delta M_{\rm disk} [M_{\odot}]$', fontsize=16)
ax1.set_title('Disk Mass - Initial Mass',fontsize=16)    

#ax2.set_ylim(0.08,11)
ax2.set_xlabel('t [yr]',fontsize=16)
ax2.set_ylabel(r'$\dot{M}_{BH} [M_{\odot}/yr]$',fontsize=16)
ax2.set_title('Black Hole Accretion Rate',fontsize=16)

ax3.set_xlabel('t [yr]',fontsize=16)
ax3.set_ylabel(r'$\dot{M}_{\rm disk} [M_{\odot}/yr]$',fontsize=16)
#ax3.set_ylim(0.0008,0.0013)
ax3.set_title('Disk Accretion Rate',fontsize=16)

plt.legend(loc='best', fontsize=16)    
plt.tight_layout()
plt.savefig('prop_all_macc.png')
plt.cla()

fig, axes = plt.subplots(figsize=(16, 12), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# choose parameter to vary here
for k in range(len(ex1)): 
    visc = vlist[1]
    mbh_init = float(mbhlist[1])*Msun
    macc_init = float(macclist[1])*Msun/yr
    ind1 = ex1[k]
    ind2 = ex2[2]
    #   General Initialization
    times = [0.0] 
    steps = 1000   
    dM_bh = [0.0]
    dM_disk = [0.0]    
    dM_acc = [macc_init]
    mdisk_init = disk.ssd_diskmass(visc, macc_init, mbh_init, ind1, ind2)
    Mbh = [mbh_init]
    Mdisk = [mdisk_init]

#   Evolution Started!
    for i in range(steps):    
    
#        macc = macc_init
#        macc = getaccrate_period(macc_init,P, amp, times[i])
        macc = getaccrate_exp(macc_init, times[i], P)
        dM_acc.append(macc)  
        a, b, c, d = diskevol(dM_acc[i+1], dM_acc[i], dM_jet, Mbh[i], Mdisk[i], dt, visc, mod, ind1, ind2)
#        b, a, d, c = disk.SSD_MEvol( Mbh[i], Mdisk[i], dM_acc[i+1], dM_acc[i], dM_jet, dt, visc, ind1, ind2)
    
        times.append( times[-1] + dt )
        dM_disk.append(a)
        dM_bh.append(b) 
        Mdisk.append(c)
        Mbh.append(d)
    
#    mod = testregion(visc,b,d)    
    # Get dt for the next step (optional)
    #dt = timestep()
       
#   Output and Plotting (for one specific case) 
       
    ax0.plot(np.array(times)/yr, (np.array(Mbh)-mbh_init)/Msun, label = 'a = %.2f'%ex1[k])
    ax1.plot(np.array(times)/yr, (np.array(Mdisk)-mdisk_init)/Msun, label = 'a = %.2f'%ex1[k])
    ax2.plot(np.array(times)/yr, np.array(dM_bh)*yr/Msun, label = 'a = %.2f'%ex1[k])    
    ax3.plot(np.array(times)/yr, np.array(dM_disk)*yr/Msun, label = 'a = %.2f'%ex1[k])
 
 
# Modification of plots etc.
ax0.set_yscale('log')
ax0.set_xlabel('t [yr]', fontsize=16)
ax0.set_ylabel(r'$\Delta M_{BH} [M_{\odot}]$', fontsize = 16)
ax0.set_title('Black Hole Mass - Initial Mass',fontsize=16)

ax1.set_yscale('log')
ax1.set_xlabel('t [yr]' ,fontsize=16)
ax1.set_ylabel(r'$\Delta M_{\rm disk} [M_{\odot}]$', fontsize=16)
ax1.set_title('Disk Mass - Initial Mass',fontsize=16)    

#ax2.set_ylim(1-0.001,1+0.001)
ax2.set_xlabel('t [yr]',fontsize=16)
ax2.set_ylabel(r'$\dot{M}_{BH} [M_{\odot}/yr]$',fontsize=16)
ax2.set_title('Black Hole Accretion Rate',fontsize=16)

ax3.set_xlabel('t [yr]',fontsize=16)
ax3.set_ylabel(r'$\dot{M}_{\rm disk} [M_{\odot}/yr]$',fontsize=16)
ax3.set_title('Disk Accretion Rate',fontsize=16)

plt.legend(loc='best', fontsize=16)    
plt.tight_layout()
plt.savefig('prop_all_a.png')
plt.cla()

fig, axes = plt.subplots(figsize=(16, 12), nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat

# choose parameter to vary here
for k in range(len(ex2)): 
    visc = vlist[1]
    mbh_init = float(mbhlist[1])*Msun
    macc_init =float(macclist[1])*Msun/yr
    ind1 = ex1[2]
    ind2 = ex2[k]
    #   General Initialization
    times = [0.0] 
    steps = 1000   
    dM_bh = [0.0]
    dM_disk = [0.0]    
    dM_acc = [macc_init]
    mdisk_init = disk.ssd_diskmass(visc, macc_init, mbh_init, ind1, ind2)
    Mbh = [mbh_init]
    Mdisk = [mdisk_init]

#   Evolution Started!
    for i in range(steps):    
    
#        macc = macc_init
#        macc = getaccrate_period(macc_init,P, amp, times[i])
        macc = getaccrate_exp(macc_init, times[i], P)
        dM_acc.append(macc)  
        a, b, c, d = diskevol(dM_acc[i+1], dM_acc[i], dM_jet, Mbh[i], Mdisk[i], dt, visc, mod, ind1, ind2)
#        b, a, d, c = disk.SSD_MEvol( Mbh[i], Mdisk[i], dM_acc[i+1], dM_acc[i], dM_jet, dt, visc, ind1, ind2)
    
        times.append( times[-1] + dt )
        dM_disk.append(a)
        dM_bh.append(b) 
        Mdisk.append(c)
        Mbh.append(d)
    
#    mod = testregion(visc,b,d)    
    # Get dt for the next step (optional)
    #dt = timestep()
       
#   Output and Plotting (for one specific case) 
       
    ax0.plot(np.array(times)/yr, (np.array(Mbh)-mbh_init)/Msun, label = 'b = %.2f'%ex2[k])
    ax1.plot(np.array(times)/yr, (np.array(Mdisk)-mdisk_init)/Msun, label = 'b = %.2f'%ex2[k])
    ax2.plot(np.array(times)/yr, np.array(dM_bh)*yr/Msun, label = 'b = %.2f'%ex2[k])    
    ax3.plot(np.array(times)/yr, np.array(dM_disk)*yr/Msun, label = 'b = %.2f'%ex2[k])
 
 
# Modification of plots etc.
ax0.set_yscale('log')
ax0.set_xlabel('t [yr]', fontsize=16)
ax0.set_ylabel(r'$\Delta M_{BH} [M_{\odot}]$', fontsize = 16)
ax0.set_title('Black Hole Mass - Initial Mass',fontsize=16)

ax1.set_yscale('log')
ax1.set_xlabel('t [yr]' ,fontsize=16)
ax1.set_ylabel(r'$\Delta M_{\rm disk} [M_{\odot}]$', fontsize=16)
ax1.set_title('Disk Mass - Initial Mass',fontsize=16)    

#ax2.set_ylim(1-0.001,1+0.001)
ax2.set_xlabel('t [yr]',fontsize=16)
ax2.set_ylabel(r'$\dot{M}_{BH} [M_{\odot}/yr]$',fontsize=16)
ax2.set_title('Black Hole Accretion Rate',fontsize=16)

ax3.set_xlabel('t [yr]',fontsize=16)
ax3.set_ylabel(r'$\dot{M}_{\rm disk} [M_{\odot}/yr]$',fontsize=16)
ax3.set_title('Disk Accretion Rate',fontsize=16)

plt.legend(loc='best', fontsize=16)    
plt.tight_layout()
plt.savefig('prop_all_b.png')
plt.cla()
