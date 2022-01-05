import numpy as np
from numba import njit
import StatsFunctions as St

@njit
def gsto_nb(rho_s, rho_o, alpha_s, alpha_o, m_a, phi, delta):
    """
        Calculate the value of the axion-photon coupling for a given value of the dark matter
        density of source and observer, the value of the stochastic parameter in source and observer,
        the mass of the ALP, the value of the amplitude on the periodic efect and a random phase. 
        Given by eq. 2.19 of arXiv 2201.xxxx 

        Parameters
        ----------
        rho_s : float
            Dark matter density at the source of the polarised light. Units of GeV/cm^3
        rho_o : float
            Dark matter density at the observer (Earth). Units of GeV/cm^3
        alpha_s : float
            Value of the alpha parameter used to parametrize time stochasticity, given by a
            Rayleigh distribution, at the source. Adimensional
        alpha_o : float
            Value of the alpha parameter used to parametrize time stochasticity, given by a
            Rayleigh distribution, at the observer. Adimensional
        m_a : float
            Mass of the axion field. Units of 10^-22 eV
        phi : float
            Value of the amplitude of the periodic effect on the polarisation of light. Units of deg
        delta : float
            Value of the delta parameter used to parametrize spatial stochasticity, given by
            a flat distribution between 0 and 2*pi. Adimensional

        Constants
        ----------
        kappa : Units conversion factor from phi (degrees) to g_a (GeV^-1)
        
        Returns
        -------
        g_agamma : float
            The coresponding axion-photon coupling for these parameters.

        Notes
        -----
        Numba @njit decorator is used to optimize its performance
    """

    kappa=6.31216*10**(-13) 
    g_agamma = kappa*m_a*phi/np.sqrt((rho_o*alpha_o**2+rho_s*alpha_s**2-2*np.sqrt(rho_o*rho_s)*alpha_o*alpha_s*np.cos(delta))/2)

    return g_agamma

@njit
def dsto_nb(rho_s, rho_o, alpha_s, alpha_o, m_a, delta_phi_aux, delta):
    """
        Calculate the value of the axion-photon coupling for a given value of the dark matter
        density of source and observer, the value of the stochastic parameter in source and observer,
        the mass of the ALP, the value of the amplitude on the periodic efect and a random phase. 
        Given by eq. 2.19 of arXiv 2201.xxxx 

        Parameters
        ----------
        rho_s : float
            Dark matter density at the source of the polarised light. Units of GeV/cm^3
        rho_o : float
            Dark matter density at the observer (Earth). Units of GeV/cm^3
        alpha_s : float
            Value of the alpha parameter used to parametrize time stochasticity, given by a
            Rayleigh distribution, at the source. Adimensional
        alpha_o : float
            Value of the alpha parameter used to parametrize time stochasticity, given by a
            Rayleigh distribution, at the observer. Adimensional
        m_a : float
            Mass of the axion field. Units of 10^-22 eV
        delta_phi_aux : float
            Value of the variance of the periodic Rayleigh distribution used to generate the
            amplitude, phi. Units of deg
        delta : float
            Value of the delta parameter used to parametrize spatial stochasticity, given by
            a flat distribution between 0 and 2*pi. Adimensional

        Constants
        ----------
        kappa : float
            Units conversion factor from phi (degrees) to g_a (GeV^-1)
        
        Returns
        -------
        delta_g_agamma : float
            The coresponding error for the axion-photon coupling for these parameters.

        Notes
        -----
        Numba @njit decorator is used to optimize its performance
    """
    kappa=6.31216*10**(-13)
    delta_g_agamma = kappa*m_a*delta_phi_aux/np.sqrt((rho_o*alpha_o**2+rho_s*alpha_s**2-2*np.sqrt(rho_o*rho_s)*alpha_o*alpha_s*np.cos(delta))/2)
    
    return delta_g_agamma

@njit
def loop(n_sources, n_universes, m_a, sigma_array, rho, n_nu):
    """
        Given the scale parameter for the Rayleigh distributions for each source and each
        frequency (mass) it generates the distribution of g_agamma resulting of the combination
        of all the sources for a number N of montecarlos (or universes).
        Explained in detail in arXiv 2201.xxxx 

        Parameters
        ----------
        n_sources : float
            Number of sources to combine on the analysis. Adimensional
        n_universes : float
            Number of Monte-Carlo simulations. Adimensiona
        m_a : numpy array
            Masses for which the analysis is done. The dimension of the array has to match
            the number of frequencies. Units of 10^-22 eV
        sigma_array : numpy array
           Scale parameter for the Rayleigh distributions for each source and each
           frequency (mass). It is an array of dimension (n_nu, n_sources). Adimensional
        rho : numpy array
            One-dimensional array of length n_sources+1 where the first entry corresponds
            to the dark matter density on earth, the rest to the dark matter density of 
            all the sources. Units of GeV/cm^3
        
        Returns
        -------
        g_global : numpy array
            Array of the value of the combination of all the sources for each MC simulation.

        Notes
        -----
        Numba @njit decorator is used to optimize its performance
    """
    g_global=[]
    for i in range(n_universes):
        
        alpha_pulsar = np.random.rayleigh(1, n_sources+1) ## stochastic case
        # alpha_pulsar = [1 for i in range(n_sources+1)] ## deterministic case  
        g_universe = []
        delta_g_universe = []
    

        for k in range(n_nu):
            
            var_pulsar =[np.complex64(x) for x in range(0)]
            g_universe_pulsars = []
            delta_g_universe_pulsars = []
            for j in range(n_sources):
                phi_pulsar = np.random.rayleigh(sigma_array[k][j])
                var_pulsar = np.sqrt(2-np.pi/2)*sigma_array[k][j]
                delta = np.random.uniform(0,2*np.pi)
                gp = gsto_nb([rho[j+1], rho[0], alpha_pulsar[j+1],alpha_pulsar[0], m_a[k], phi_pulsar,delta])
                dp = dsto_nb([rho[j+1], rho[0], alpha_pulsar[j+1],alpha_pulsar[0], m_a[k], var_pulsar,delta])
                g_universe_pulsars.append(gp)
                delta_g_universe_pulsars.append(dp)

                
            mw = St.meanw_nb(g_universe_pulsars,delta_g_universe_pulsars)[0]    
            uw = St.meanw_nb(g_universe_pulsars,delta_g_universe_pulsars)[1]    
            g_universe.append(mw)
            delta_g_universe.append(uw)
        g_global.append(g_universe)
    return g_global


def NuToMa(nu):
    ma = 1.31/(2*np.pi / nu / 365)
    return  ma