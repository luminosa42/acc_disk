# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import disk

G = 6.67259 * 10.**-8 
cc = 2.99792458 * 10.**10
Msun = 1.9891 * 10.**33  
yr = 3600. * 24. * 365.25  
pc = 3.086 * 10. ** 18
e = 0.1

def get_tacc(eta, mbh, visc, a, b):
    
    K = disk.ssd_diskconst(visc)
    cedd = 1.25e38/ cc **2./ Msun /e
    tacc = K * a * cedd ** (a-1) * mbh ** (a+b-1) * eta **(a-1)
    return tacc
    
def get_cs(mbh, rdisk, hr):    

    return hr* (G * mbh/ rdisk) **(1./2.)
    
def get_sound_damp(mbh, rdisk, hr):

    cs = get_cs(mbh, rdisk, hr)
    tsound = 5. * 10**5. * yr * rdisk/(0.1* pc) * 2.* 10**4. / cs
    return tsound
    
def get_tvisc(visc, mbh, rdisk, hr):
    
    tvisc = 2./3./visc * hr **(-2.) * (G * mbh/ rdisk**3.) **(-1./2.)
    return tvisc

# inputs
vlist = [0.01, 0.1, 0.25]
mbhlist = [1.e8, 1.e9, 1.e10]
ex1 = [1./27., 1./9., 5./27.]
ex2 = [1./3., 1./2., 2./3.]    
plt.figure(figsize=(12,9)) 

for i in range(len(mbhlist)):
    
    mbh = float(mbhlist[i])*Msun
    visc = vlist[1]
    a = ex1[2]
    b = ex2[2]
    etarange = np.linspace(-1,1,100)    # this is in log scale
    tacclist = []
    tsoundlist = []
    tvisclist = []
    
    for eta in etarange:
    
        eta = 10**eta
        tacclist.append(get_tacc(eta,mbh,visc,a,b))
        macc = eta * disk.Medd(mbh)
        rdisk = disk.ssd_Rmax(macc, mbh, visc)
        hr = disk.ssd_hr(macc, mbh, rdisk, visc)
        print (macc*yr/Msun, rdisk/pc, hr)
        tsoundlist.append(get_sound_damp(mbh,rdisk,hr))
        tvisclist.append(get_tvisc(visc,mbh,rdisk,hr))

    plt.plot(etarange, np.array(tacclist)/yr, label = r'$t_{\rm acc}\ M_{\rm BH} = %e M_{\odot}$'%mbhlist[i])
    plt.plot(etarange, np.array(tsoundlist)/yr, linestyle='--',label = r'$t_s\ M_{\rm BH} = %e M_{\odot}$'%mbhlist[i])
    plt.plot(etarange, np.array(tvisclist)/yr, linestyle='-.', label = r'$t_{\rm visc}\ M_{\rm BH} = %e M_{\odot}$'%mbhlist[i])
   
   
plt.yscale('log')
plt.xlabel(r'$\log(\dot{M}/\dot{M}_{edd})$', fontsize = 16)
plt.ylabel(r'$\rm t [yr] $', fontsize = 16)
#plt.ylim(1,1e8)
plt.legend(loc='best', fontsize = 16)
plt.title(r'$\alpha = 0.1$', fontsize = 20)
plt.savefig('tacc.png')    