# -*- coding: utf-8 -*-

import numpy as np

Msun = 1.9891 * 10.**33
G = 6.67259 * 10.**-8                 
cc = 2.99792458 * 10.**10  
e = 0.1
#e = 1.
beta = 0.5
#-------------------------------------------------------------------------------
#---------------------- GENERAL FUNCTIONS FOR ALL DISKS ------------------------
#-------------------------------------------------------------------------------

def Ledd(M):   # Eddington Luminosity
    
    return ( 1.3 * 10.** 38) * M/Msun
    
def Rsch(M): # Schwarzchild radius

    return 2.0 * G * M / cc**2.
    
def Risco(M) : # Innermost stable circular orbit

    return 6.0 * G * M / cc**2.   
    
def Medd(M):
    
    return Ledd(M)/(e * cc**2.)

#-------------------------------------------------------------------------------
#------------------------ Place for any test functions --------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#----------------------- SHAKURA SUNYAEV THIN DISK -----------------------------
#-------------------------------------------------------------------------------
def ssd_diskconst(visc):
    
    return (3.58E6) * (visc/0.03)**(-2./27.)    
    
def ssd_diskmass(visc, Macc, Mbh, a, b):
    
    const = ssd_diskconst(visc)
    Mdisk = const * ( Macc**a ) * ( Mbh**b)
    return Mdisk

def ssd_diskacc(visc, macc_new, macc_old, dt, Mbh, dM_jet, a, b):
    
    const = ssd_diskconst(visc)
    timedep = const * a * (macc_new **(a-1.) ) * (Mbh ** b) * (macc_new-macc_old)/dt
    dM_disk = const * b * ( macc_new**(a+1.) ) * ( Mbh**(b-1.) ) - const * b * ( macc_new**a ) * ( Mbh**(b-1.) ) * dM_jet + timedep
    dM_disk = dM_disk/ (1. + const * b * ( macc_new**a ) * ( Mbh**(b-1.) ) )
    return dM_disk

def ssd_hr(Macc, Mbh, R, visc):
    
    #Calculates and returns H/R at radius R; unitless.
    Rs = Rsch(Mbh)
    Le = Ledd(Mbh)
    L = e * Macc * cc**2
#    R = ssd_Rmax(Macc, Mbh, visc)

    return ( 1.94*10.**-3 ) * ( visc/0.03 )**(-0.1) * ( e/0.1 )**(-0.2) * \
           ( L/(0.1*Le) )**(0.2) * (Mbh/10.**8 /Msun )**(-0.1) * ( R/Rs )**(0.05)

def ssd_Rmax(Macc, Mbh, visc):
    
    Le = Ledd(Mbh)
    L = e * Macc * cc**2
    Rs = Rsch(Mbh)
    #Calculates and returns outer radius of the disk a as R/Rs
    return ( 1.13*10.**3 )  * ( visc/0.03 )**(14./27.) * ( e/0.1 )**(8./27.) * \
           ( L/(0.1*Le) )**(-8./27.) * (Mbh/10.**8 /Msun )**(-26./27.) * Rs
           
def ssd_surfdens(Macc, Mbh, R, visc ):
    
    Le = Ledd(Mbh)
    L = e * Macc * cc**2
    Rs = Rsch(Mbh)
    #Calculates and returns surface density at R; units of g cm^-2
    return ( 7.5*10.**6 ) * ( visc/0.03 )**(-0.8) * ( e/0.1 )**(-0.6) *   \
           ( L/(0.1*Le) )**(0.6) * (Mbh/10.**8 /Msun )**(0.2) * ( R/Rs )**(-0.6) 
           
#-------------------------------------------------------------------------------
#----------------- Optically-thick, two-temperature ADAF -----------------------
#-------------------------------------------------------------------------------
def adaf_diskconst(visc):

    return 1.21E6 * 4 * np.pi *2./3. /visc
    
def adaf_rmax(Macc, Mbh, visc):
    
    Le = Ledd(Mbh)
    return ( 2.1*10.**3 ) * ( visc/0.3 * beta * Le**2. * (e/0.1) **-2. / (Mbh/10.**8 /Msun ))**(2./9.)* Rsch(Mbh)
           
def adaf_diskmass(visc, Macc, Mbh, a, b):
    
    const = adaf_diskconst(visc)
    Me = Medd(Mbh)
    return const * (Mbh/Msun) ** -1. * (Macc/Me) 
    
def adaf_diskacc():
    
    return 0
    
    
    
    
    