"""
Code to produce multiple HR diagrams given RGB fluxes, distance compared with baseline. Compare with baselineHR.ipynb
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import os 

def fluxToApp(flux):
    return np.log(flux)

def appToAbs(m,d):
    # Convert apparent to absolute magnitudes given distance.
    return m + 2*np.log(d) 

def plotHR(MRs, MGs, MBs, MPs, ax: plt.Axes, col=None, lab='', scale=1):
    """
    Given a collection (iterables of same length) of red, green, blue
        absolute magnitudes, compute & plot the HR diagram, 
        on the given axis "ax" with colour "col" and label "lab".
    Note that this does not create figures.
    """
    # plot all stars with no listed period
    class0Stars = MPs == 0
    ax.scatter(MRs[class0Stars] - MBs[class0Stars], MGs[class0Stars], marker='.', 
        c=col, s=scale, label=lab)

    # plot 18-21 period (class 1) stars with diamond
    # take bool arrays, cast to int, then multiply
    class1Stars = (18 < MPs).astype('int') * (MPs < 21).astype('int') 
    ax.scatter(MRs[class1Stars] - MBs[class1Stars], MGs[class1Stars], marker='x', 
        c=col, s=scale, label=lab + ' Class 1')

    # plot 48-52 period (class 2) stars with square
    class2Stars = (48 < MPs).astype('int') * (MPs < 52).astype('int') # bool array of stars class 2
    ax.scatter(MRs[class2Stars] - MBs[class2Stars], MGs[class2Stars], marker='^', 
        c=col, s=scale, label=lab + ' Class 2')

    # there SHOULD NOT be other stars that do not fit in these three classes
    ax.grid()

def plotBaseline(colour='b'):
    """
    c.f. baselineHR.ipynb for explanation, this is here
        for convenience in later plotting
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    calStars = pd.read_csv(dir_path + '/../Data/Camera Images/calibrated star data.txt',delimiter=' ')
    
    parallaxLimit = 0.01

    filtStars = calStars.loc[calStars['Parallax'] > parallaxLimit]
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
            [filtStarRed, filtStarGreen, filtStarBlue])), filtStarPer, ax, col=colour, lab='Benchmark')

    ax.set_title('HR Diagram')
    ax.set_ylabel('Absolute green magnitudes')
    ax.set_xlabel('$\log(F_{Red}) - \log(F_{Blue})$')

    return fig, ax

if __name__ == "__main__":
    # put code that you want to run from this file here
    warnings.filterwarnings('ignore')

    # current code is just to make sure it works properly!
    fig, ax = plotBaseline()
    fig.show()
