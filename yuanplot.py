# -*- coding: utf-8 -*-

import numpy as np
import disk
import main2
import matplotlib.pyplot as plt
import scipy.optimize as op

# Physical Constants
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
gamma = 5./3.
f = 1.
sadaf = 0.
N = 200

m = 10.
Mbh = m * Msun
macc_array = np.linspace(-5, 5, N)
visc = 0.01
    
mod = []
dens = []
Tmp = []
hline = [] # where h/r = 1
alfsig1 = []
alfsig2 = []

Nr = 5.
r = Nr * disk.Rsch(Mbh)

def hr_eqns(x, consta, constb ):
    
    [rho,temp] = x
    pres = main2.pressure(rho, temp)
    cs = np.sqrt(pres/rho)
    
    return [cs- consta,]

for i in range(N): 
    
    macc_scale = 10. ** macc_array[i]
    mbh_scale = m
    Macc = macc_scale * disk.Medd(Mbh)
    dr = disk.Risco(Mbh)
    vk = main2.orbitalspeed(r, Mbh)
     
    # the hr curve
    
#    hline.append(np.log10(2. * dens[i]* r* visc))
    modi = main2.testregion(visc, Nr, mbh_scale, macc_scale)
#    modi = 'ssd'
    mod.append(modi)
    
    if (mod[i] == 'ssd'):
                
        [rho, tmp, pi] = main2.ssd_regions(visc, Nr, dr, Mbh, Macc)
        
        dens.append(rho)
        Tmp.append(tmp)
        
        #pi = main2.pressure(dens[i], Tmp[i])   
        #pi = dens[i]* kb * Tmp[i]/ mu/ mH                           
        csi = np.sqrt(pi/dens[i])        
        hi = csi/vk
        alfsig1.append(np.log10(2. * dens[i]* hi* visc))
        alfsig2.append(np.nan)
        
        
    elif(mod[i] == 'adaf'):
        
        """
        epsilon = (5./3. - gamma)/(gamma-1.)/f
        g = np.sqrt(1. + 18. * visc **2. /(5. + 2. * epsilon) **2. ) - 1.
        vff = main2.ffspeed(r, Mbh)
        vri = (5. + 2. * epsilon)/(3. * visc) * g * vff
        csi = np.sqrt(2. * (5. + 2. * epsilon)/(9. * visc ** 2.) * g * vff **2.)
        hi = np.sqrt(2.5) * csi/vff * r     
        rhoi = Macc/(4. * np.pi * r * hi * vri)
        dens.append(rhoi)
        pi = rhoi * csi **2.
        """
        
        rscale = Nr
        csi = np.sqrt(1.4e20/rscale)
        vri = 1.1e10 * visc * rscale ** (-1./2.)
        pi = 1.7e16 * macc_scale / (visc * mbh_scale) * rscale ** (-5./2.+sadaf)
        omegai = 2.9e4/ (mbh_scale * rscale ** (3./2.))
        rhoi = pi/ csi **2. 
        hi = csi/omegai       
        dens.append(rhoi)

        adaf_sol = op.root(main2.adaf_temp, 1e4, args = (rhoi, pi,))        
        Tmp.append(adaf_sol.x)
        
        alfsig2.append(np.log10(2. * dens[i]* hi* visc))
        alfsig1.append(np.nan)
        
    else:
        
        Tmp.append(0)
        dens.append(0)
        alfsig1.append(np.nan)
        alfsig2.append(np.nan)

#print(alfsig1)        

# the M-sigma plot     

plt.plot(alfsig1, macc_array, label = 'SSD', linewidth=2)
plt.plot(alfsig2, macc_array, label = 'ADAF', linewidth=2)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.axvline(x=0., color='k')
plt.xlabel(r'$\log (\alpha \Sigma)$')
plt.ylabel(r'$\log(\dot{M}/\dot{M}_{\rm Edd})$')
plt.title(r'$M/M_{\odot} = %.1e, r = 5$'%m)
plt.legend()
plt.savefig('equl.png')
plt.cla()

# the h-R plot

    
