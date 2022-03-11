# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 11:02:06 2022

@author: ryanw
"""

from numpy import *
from astropy.timeseries import LombScargle
from astropy.table import Table
import matplotlib.pyplot as plt
import os 


dir_path = os.path.dirname(os.path.realpath(__file__))  #finds the path of this program to use later

stars = Table.read(dir_path + "/calibrated star data.txt", format='ascii')

mb, mg, mr = (log10(stars['BlueFlux']), log10(stars['GreenFlux']), log10(stars['RedFlux']))
colour = mr - mb
absmag = mg + 2 * log10(stars['Distance'])

#ax1 = plt.subplots(1,1)

plt.scatter(colour, absmag)