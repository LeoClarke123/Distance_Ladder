# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:43:01 2022

@author: Ryan White     s4499039

Plots a basic star map of the sky as seen from the observer. 
There's quite a bit i'm not happy with in here, but sunk cost fallacy made me get it to a useable state. 

Needs fixing:
    - The brightness of each point
    - The colour of each point (they're currently too vibrant)
    - Each scatter point has a strange circle around it which I'm not sure how to fix.
"""

from numpy import *
from astropy.timeseries import LombScargle
from astropy.table import Table
import os 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


dir_path = os.path.dirname(os.path.realpath(__file__))  #finds the path of this program to use later
totalpoints = open(dir_path+"/"+"total point-like data.txt", "r")       #opens and reads existing data. 

data = totalpoints.readlines()[1:] 
allaveflux = []
stardata = []
allequats, allpolar, allbright, colours = [], [], [], []

for row in data:
    [name, equat, polar, bluef, greenf, redf, parallax, veloc, distance, period, location] = row.split(", ")
    name, location = name.replace('[', ''), location.replace(']', '')
    bluef, greenf, redf = float(bluef), float(greenf), float(redf)
    rgb = [bluef, greenf, redf]
    aveflux = average(rgb)
    allaveflux.append(aveflux)
    stardata.append([equat, polar, aveflux])
    allequats.append(float(equat))
    allpolar.append(180 - float(polar))         #flips the polar coords
    norm = max(rgb)         #finds max flux value for normalization purposes.
    r, g, b = 255 * redf / norm, 255 * greenf / norm, 255 * bluef / norm
    colours.append(array([r, g, b]))
    
brightest = max(allaveflux)         #finds the brightest star for normalization purposes

for star in stardata:
    brightness = sqrt(11) * log10((star[2] / brightest) + 1)            #defines scatter marker 'brightness' based on logarithmic algorithm. Spits out value between 0 and 1
    allbright.append(brightness)

colours = array(colours)/256        #gets the rgb values between 0 and 1
plt.scatter(allequats, allpolar, s=allbright, c=colours, marker='.')
ax = plt.gca()          #i don't even know what this does
fig = plt.gcf()
plt.xlabel('Equatorial Angle (deg)')
plt.ylabel('Inverse Polar Angle (deg)')
ax.set_facecolor('k')
figure(figsize=(36,18))             #units are inches
fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
fig.savefig('starmap.png')
