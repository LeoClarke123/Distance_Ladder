# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:11:16 2022

@author: Ryan White     s4499039

This function takes star spectra data for stars of known distance, and attempts to plot a planck curve to the data. 
Once this has been done, the temperature of the star is taken and then all of the temperatures of the stars are plotted on a histogram. 
"""
import statistics
from numpy import *
import os 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import logging

logging.basicConfig(level=logging.INFO)
matplotlib.use('Agg') # ! (Ciaran)  can't run the code without this, not sure why


def planck_func(x, T, a):
    '''This is a planck function to be used with the scipy.optimize.curve_fit() function. 
    '''
    const1 = 4.702 * 10**(-15)         #corresponds to 8 pi^2 * h * c^2
    const2 = 0.01439              #corresponds to h * c / k 
    funct = a * (const1 / (x)**5) * (1 / (exp(const2 / (x * T)) - 1))
    return funct

dir_path = os.path.dirname(os.path.realpath(__file__))

stardata = pd.read_csv(dir_path + "\\Camera Images\\calibrated star data.txt", delimiter=' ')
starLabel = stardata['Name']
starRed = stardata['RedFlux']; starGreen = stardata['GreenFlux']; starBlue = stardata['BlueFlux']
starDist = stardata['Distance']

logging.debug(f"len(starRed) = {len(starRed)}")

r, g, b = 700 * 10**(-9), 500 * 10**(-9), 440 * 10**(-9)        #corresponding wavelengths for each luminosity measurement
wavelengths = array([r, g, b], dtype=float64)
colours = ["r", "g", "b"]

temps = []; lowtemps = []; names = []
uncertainties = []
for index, star in enumerate(starRed):
    distMult = (starDist[index] * 3.086 * 10**16)**2        #approximates the actual luminosity by multiplying the lumin by the distance to it
    lumins = array([starRed[index]*distMult, starGreen[index]*distMult, starBlue[index]*distMult], dtype=float64)   #RGB luminosity for the star
    try:
        fig, ax = plt.subplots()
        
        starProps, params = curve_fit(planck_func, wavelengths, lumins, p0=[10000, 5 * 10**7])      #fits the planck curve to the data, with initial guess T=10000K
        starTemp, mult = starProps
        if starTemp < 200000:
            names.append(starLabel[index])
            temps.append(starTemp)
            uncertainties.append(sqrt(diag(params))[0])
        if starTemp < 50000:        #for the zoomed in histogram
            lowtemps.append(starTemp)
        # if sqrt(diag(params))[0] < 200:
        x = arange(10*10**(-9), 1300*10**(-9), 900*10**(-9) / 50)
        planck = mult * (4.702 * 10**(-15) / (x)**5) * (1 / (exp(0.01439 / (x * starTemp)) - 1))
        plt.plot(x * 10**9, planck, c='k', linewidth=0.7)
        for i in range(0,3):
            plt.scatter(wavelengths[i] * 10**9, lumins[i], c=colours[i], s=10)
            plt.errorbar(wavelengths[i] * 10**9, lumins[i], yerr=0.015*lumins[i], linewidth=1, fmt=',', c=colours[i])
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Specific Luminosity (W.nm$^{-1}$)")
        
        ax.ticklabel_format(axis='x', style='scientific', useMathText=True)
        fig.savefig(dir_path+f'\\Blackbody Curves\\{stardata["Name"][index]}-T{round(starTemp)}K.png', dpi=100, bbox_inches='tight', pad_inches = 0.01)
        fig.savefig(dir_path+f'\\Blackbody Curves\\{stardata["Name"][index]}-T{round(starTemp)}K.pdf', dpi=100, bbox_inches='tight', pad_inches = 0.01)
        plt.close(fig)
    except RuntimeError:        #for cases where the curve_fit function can't find the temperature
        pass
    finally:
        plt.close(fig)

# write star label, temp, uncertainty to a csv
pd.DataFrame([names, temps, uncertainties]).transpose().to_csv(dir_path+'\\temperatureUncertainties.csv')

#this plots all of the startemps in a histogram
fig, ax = plt.subplots()
plt.hist(temps, bins=80)
ax.set_xlabel("Star Temperature (K)")
ax.set_ylabel("Number of Nearby Stars")
fig.savefig(dir_path+'\\AllStarTemps.png', dpi=400, bbox_inches='tight', pad_inches = 0.01)
fig.savefig(dir_path+'\\AllStarTemps.pdf', dpi=400, bbox_inches='tight', pad_inches = 0.01)
plt.close(fig)


#the following plots the main-sequence range star temps in a histogram with a cumulative distribution overlay
fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(lowtemps, bins=50)
ax.grid(which='major', axis='y')
ax.set_xlabel("Star Temperature (K)")
ax.set_ylabel("Number of Nearby Stars")

values, base = histogram(lowtemps, bins=50)
cumulative = cumsum(values)
ax2 = ax.twinx()
ax2.plot(base[:-1], cumulative / max(cumulative), c='r')
plt.ylim([0,1.05])
ax2.set_ylabel("Cumulative Distribution of Star Temp", c='maroon')

fig.savefig(dir_path+'\\MainSeqStarTemps.png', dpi=400, bbox_inches='tight', pad_inches = 0.01)
fig.savefig(dir_path+'\\MainSeqStarTemps.pdf', dpi=400, bbox_inches='tight', pad_inches = 0.01)
plt.close(fig)

logging.info(f"Low temps median {statistics.median(lowtemps)}")
