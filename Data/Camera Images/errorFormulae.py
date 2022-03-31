"""
File for analysing the uncertianty in variable star periods using Monte Carlo estimation

Add your uncertainty formulae in here!

@author: leotc & ryanw
"""

import numpy as np
from astropy.timeseries import LombScargle
from astropy.table import Table

def period_uncert(star_time, star_flux, runs = 1000,  yerr = 0.015):
    """
    Function for calcuationg the uncertianty in periods of variable stars using monte carlo method
    
    Inputs: star_time <astropy.table>: variable star recording times
            star_flux <astropy.table>: variable star flux readings
            runs <int>: number of monte carlo iterations to perform, defaul = 1000
            yerr <float>: percentage uncertainty in flux readings, default = 1.5%            
    
    Returns: <float>: uncertainty in variable star period (standard deviation of simulation)
    """
    np.random.seed(1001)
    
    i=0
    var_period = []
    
    while i<runs:
        
        # add random error to each flux term based on uniform distribution of uncertainties
        flux_jitter = star_flux + (star_flux*yerr*np.random.uniform(low=-1, high=1, size=len(star_flux)))

        LS = LombScargle(star_time, flux_jitter)
        freqs = np.linspace(1/100, 0.45, 10000)        #frequency grid
        power = LS.power(freqs)     #calculates LS power
        var_period.append(1 / freqs[np.argmax(power)])
        
        i+=1
    
    # return standard deviation of simulation as uncertianty
    return np.std(var_period) 