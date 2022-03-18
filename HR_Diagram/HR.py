"""
Code to produce multiple HR diagrams given RGB fluxes, distance compared with baseline. Compare with baselineHR.ipynb
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import os 


def appToAbs(m,d):
    # Convert apparent to absolute magnitudes given distance.
    return m + 2*np.log(d) 

def plotHR(MRs, MGs, MBs, ax: plt.Axes, col=None):
    """
    Given a collection (iterables of same length) of red, green, blue
        absolute magnitudes, compute & plot the HR diagram, 
        on the given axis "ax" with colour "col".
    Note that this does not create figures.
    """
    ax.scatter(MRs - MBs, MGs, edgecolors=col)

    ax.set_title('HR Diagram')
    ax.set_ylabel('Absolute green magnitudes')
    ax.set_xlabel('$\log(F_{Red}) - \log(F_{Blue})$')
    ax.grid()

def plotBaseline():
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

    fig, ax = plt.subplots(1,1)
    plotHR(*list(map(lambda x: appToAbs(np.log(x), filtStarDists),[filtStarRed, filtStarGreen, filtStarBlue])), ax)

    return fig, ax

if __name__ == "__main__":
    # put code that you want to run from this file here
    warnings.filterwarnings('ignore')

    # current code is just to make sure it works properly!
    fig, ax = plotBaseline()
    fig.show()
