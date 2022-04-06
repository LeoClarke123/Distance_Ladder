"""
Code to produce multiple HR diagrams given RGB fluxes, 
distance compared with baseline. Compare with baselineHR.ipynb.
Confer distances.pdf for the units used in these diagrams.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import os

__author__ = "Ciaran Komarakul-Greene, s4528583"

PARALLAX_LIMIT = 0.01
UNCERTAINTY_LOWER = 0.78
UNCERTAINTY_UPPER = 1.29

def fluxToApp(flux):
    return -np.log(flux)

def appToAbs(m,d):
    # Convert apparent to absolute magnitudes given distance.
    return m - 2*np.log(d)

def distModToDist(m_M):
    return np.exp((m_M)/2)
    
def distShiftByModulus(d,m):
    # shift a distance in pc given a modulus to shift UP by
    # to be used to make adjustments to HR diagrams
    return distModToDist(2*np.log(d)+m)

def plotHR(MRs, MGs, MBs, MPs, ax: plt.Axes, col=None, lab='Cluster ', scale=1):
    """
    Given a collection (iterables of same length) of red, green, blue
        absolute magnitudes, compute & plot the HR diagram, 
        on the given axis "ax" with colour "col" and label "lab".
    Note that this does not create figures.
    """
    # plot all stars with no listed period
    class0Stars = MPs == 0
    ax.scatter(MBs[class0Stars] - MRs[class0Stars], MGs[class0Stars], marker='.', 
        c=col, s=scale, label=lab)

    # Note that we double the scale of variable stars to make them more pronounced
    # Plot 18-21 period (class 1) stars with diamond marker.
    class1Stars = np.logical_and(18 < MPs, MPs < 21) 
    ax.scatter(MBs[class1Stars] - MRs[class1Stars], MGs[class1Stars], marker='x', 
        c=col, s=2*scale, label=lab + ' Class 1')

    # plot 48-52 period (class 2) stars with square marker
    class2Stars = np.logical_and(48 < MPs, MPs < 52) 
    ax.scatter(MBs[class2Stars] - MRs[class2Stars], MGs[class2Stars], marker='^', 
        c=col, s=2*scale, label=lab + ' Class 2')

    # there SHOULD NOT be other stars that do not fit in these three classes
    ax.grid()

def plotBaseline(colour='b'):
    """
    c.f. baselineHR.ipynb for explanation, this is here
        for convenience in later plots
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    calStars = pd.read_csv(dir_path + '/../Data/Camera Images/calibrated star data.txt',
        delimiter=' ')

    filtStars = calStars.loc[calStars['Parallax'] > PARALLAX_LIMIT]
    filtStarBlue = np.array(filtStars['BlueFlux']) # blue flux
    filtStarGreen = np.array(filtStars['GreenFlux']) # green flux
    filtStarRed = np.array(filtStars['RedFlux']) # red flux
    filtStarDists = np.array(filtStars['Distance']) # distances
    filtStarPer = np.array(filtStars['Periodicity']) # periods

    fig, ax = plt.subplots(1,1)
    # the first argument to plotHR does the following:
    #   maps appToAbs on each array of fluxes given the distance array filtStarDists
    #   then converts to list, and unpacks as argument to plotHR 
    plotHR(*list(map(lambda x: appToAbs(fluxToApp(x), filtStarDists),
            [filtStarRed, filtStarGreen, filtStarBlue])), filtStarPer, ax, col=colour, 
            lab='Benchmark')

    ax.set_ylabel('$M_G$')
    ax.set_xlabel('$M_B - M_R$')
    ax.invert_yaxis()

    return fig, ax

def main():
    print("This code doesn't do anything by itself, why are you running it")
    # I am funny please laugh

if __name__ == "__main__":
    # put code that you want to run from this file here
    warnings.filterwarnings('ignore')

    main()
