# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 10:45:01 2022

@author: leotc
"""

from numpy import *
from math import isclose
import os 
import statistics
import matplotlib as mpl


dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the Cluster Identifier.py directory
totalpoints = open(dir_path + "/calibrated star data.txt", "r")      
pointdata = totalpoints.readlines()[1:]

distance = []
period = []
luminosity = []

for i in pointdata:
    i = i.split()
    if float(i[-2]) != 0:
        period.append(float(i[-2]))
        abs_mag = np.log10(float(i[4])) + 2*np.log10(float(i[-3]))
        luminosity.append(abs_mag)


mpl.pyplot.scatter(period, luminosity, marker='.')
mpl.pyplot.xlabel('Period (h)')
mpl.pyplot.ylabel('Log Flux')