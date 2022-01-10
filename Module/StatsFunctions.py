import numpy as np
from numba import njit

def sigray(phi_95):
    """
        Calculate the scale parameter sigma of a rayleigh distribution
        from the 95% value of its CDF.

        Parameters
        ----------
        phi_95 : float
            The value of the CDF at 95% confidence interval for the rayleigh distribution.

        Returns
        -------
        sr : float
            The scale parameter, sigma, of a rayleigh distribution.

        Notes
        -----
        The 95% confidence interval for the standard deviation of a rayleigh
        distribution is given by:
        
        .. math::
            \\sqrt{\\frac{\\phi^2}{\\left|2\\log(0.05)\\right|}}

        where :math:`\\phi` is the 95% confidence interval for the rayleigh
        distribution.
    """
    sr = np.sqrt(phi_95**2/np.abs(2*np.log(0.05)))
    return sr

def medray(phi_95):##median
    """
        Calculate the median of a rayleigh distribution
        from the 95% value of its CDF.

        Parameters
        ----------
        phi_95 : float
            The value of the CDF at 95% confidence interval for the rayleigh distribution.

        Returns
        -------
        mr : float
            The median of a rayleigh distribution.

        Notes
        -----
        The median of a rayleigh distribution is given by:
        
        .. math::
            \\sigma\\sqrt{2\\ln(2)}

        where :math:`\\sigma` is the scale parameter for the rayleigh
        distribution.
    """
    mr = np.sqrt(phi_95**2/np.abs(2*np.log(0.05)))*np.sqrt(2*np.log(2))
    return mr

@njit
def meanw_nb(values, weights):
    """ Calculate the weighted mean of a list of values.
        Parameters
        ----------
        values : list of float
            The list of values to be averaged.
        weights : list of float
            The list of weights for each value.

        Returns
        -------
        mean : float
            The weighted mean of the input values.
        uncertainty : float
            The uncertainty in the weighted mean."""
    aw = np.sum(np.array(values)/np.array(weights)**2)
    va = np.sum(1/np.array(weights)**2)
    mw = aw / va
    return mw, np.sqrt(1/va)
