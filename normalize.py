import numpy as np

def getSelfSimilarValues (mvir, delta, crit=True, aexp=1.0, omega_m=0.27, omega_l=0.73,	omega_b = 0.0469, hubble=0.7 ) :
    ''' 
    Return self-similar quantites (e.g. T500c) for temperature [keV] pressure [erg/cm^3], and entropy [keV cm^2] given halo mass mvir [Msun/h], its overdensity delta, with respect to critical (crit=True) or mean (crit=False) at expansion factor aexp (= 1.0 by default). Assume WMAP5 cosmological parameters (true for Bolshoi and L500 simulations).
    '''
    fb = omega_b/omega_m
    Ez = np.sqrt(omega_m/aexp**3.0+omega_l)
    mu = 0.59
    mue = 1.14

    mpivot = 1e15 # in Msun/h

    delta = float(delta)
	
    if crit :
        Tfit_norm = 11.055*(mu/0.59)*(hubble/0.7)**(2./3.)
        Tfit_alpha = 2./3.
        Tfit_beta = 2./3.
        Pfit_norm = 1.4458e-11*(fb/0.1737)*(hubble/0.7)**2
        Pfit_alpha = 2./3.
        Pfit_beta = 8./3.
        Kfit_alpha = 2./3.
        Kfit_beta = -2./3.
        #Kfit_norm = 1963.6*(mu/0.59)*(mue/1.14)**(2./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)
        Kfit_norm = 1265.7*(mu/0.59)**(5./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)

        Tvir = (delta/500.0)**(1./3.)*Tfit_norm*(mvir/mpivot)**(Tfit_alpha)*Ez**(Tfit_beta)
        Kvir = (delta/500.0)**(-1./3.)*Kfit_norm*(mvir/mpivot)**(Kfit_alpha)*Ez**(Kfit_beta) 
        Pvir = (delta/500.0)**(4./3.)*Pfit_norm*(mvir/mpivot)**(Pfit_alpha)*Ez**(Pfit_beta)
    else :
        Tfit_norm = 5.2647*(mu/0.59)*(hubble/0.7)**(2./3.)*(omega_m/0.27)**(1./3.)
        Tfit_alpha = 2./3.
        Tfit_beta = -1.0
        Pfit_norm = 7.4359e-13*(fb/0.1737)*(hubble/0.7)**2*(omega_m/0.27)**(4./3.)
        Pfit_alpha = 2./3.
        Pfit_beta = 4.0
        Kfit_alpha = 2./3.
        Kfit_beta = 1.0
        #Kfit_norm = 4123.3*(mu/0.59)*(mue/1.14)**(2./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)*(omega_m/0.27)**(-1./3.)
        Kfit_norm = 2799.75*(mu/0.59)**(5./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)*(omega_m/0.27)**(-1./3.)
        Tvir = (delta/200.0)**(1./3.)*Tfit_norm*(mvir/mpivot)**(Tfit_alpha)*(aexp)**(Tfit_beta)
        Kvir = (delta/200.0)**(-1./3.)*Kfit_norm*(mvir/mpivot)**(Kfit_alpha)*(aexp)**(Kfit_beta)
        Pvir = (delta/200.0)**(4./3.)*Pfit_norm*(mvir/mpivot)**(Pfit_alpha)/(aexp)**(Pfit_beta)

    return (Tvir, Pvir, Kvir)