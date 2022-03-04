# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 14:43:01 2022

@author: Ryan White     s4499039

Plots a basic star map of the sky as seen from the observer. Also produces a basic fuzzy-object map in the same format.
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
totalpoints = open(dir_path+"/total point-like data.txt", "r")       #opens and reads existing data. 
totalfuzzy = open(dir_path+"/total fuzzy data.txt", "r")


stardata = totalpoints.readlines()[1:]      #assigns the data
fuzzydata = totalfuzzy.readlines()[1:]

allstaraveflux, allfuzaveflux = [], []      #initialize arrays
stararray, fuzarray = [], []
starequats, starpolar, starbright, starcolours = [], [], [], []
fuzequats, fuzpolar, fuzbright, fuzcolours = [], [], [], []

for row in stardata:
    [name, equat, polar, bluef, greenf, redf, parallax, veloc, distance, period, location] = row.split()
    bluef, greenf, redf = float(bluef), float(greenf), float(redf)       #cleans up variables
    rgb = [bluef, greenf, redf]
    aveflux = average(rgb)          #averages the r g b flux values
    allstaraveflux.append(aveflux)      #adds it to a list
    stararray.append([equat, polar, aveflux])       
    starequats.append(float(equat))
    starpolar.append(float(polar))         #flips the polar coords
    norm = max(rgb)         #finds max flux value for normalization purposes.
    r, g, b = 255 * redf / norm, 255 * greenf / norm, 255 * bluef / norm            #assigns the rgb values
    starcolours.append(array([r, g, b]))
    
brightest = max(allstaraveflux)         #finds the brightest star for normalization purposes

for star in stararray:
    brightness = sqrt(11) * log10((star[2] / brightest) + 1)            #defines scatter marker 'brightness' based on logarithmic algorithm. Spits out value between 0 and 1
    starbright.append(brightness)

starcolours = array(starcolours)/256        #gets the rgb values between 0 and 1 (because the scatter function chucked a tantrum)

fig, ax = plt.subplots()            #initialize axes
plt.xlabel('Equatorial Angle (deg)')
plt.ylabel('Polar Angle (deg)')
ax.set_facecolor('k')
plt.gca().invert_yaxis()        #flips the y-axis so that polar angle is 0 at due north

plt.scatter(starequats, starpolar, s=starbright, c=starcolours, marker='.')

figure(figsize=(36,18))             #units are inches
fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
fig.savefig('starmap.png')

totalpoints.close()
plt.clf()           #clears the current figure in order to create the next one

for row in fuzzydata:       #functionally identical to the star loop
    [name, equat, polar, bluef, greenf, redf, size, veloc, location] = row.split()
    bluef, greenf, redf = float(bluef), float(greenf), float(redf)
    rgb = [bluef, greenf, redf]
    aveflux = average(rgb)
    allfuzaveflux.append(aveflux)
    fuzarray.append([equat, polar, aveflux])
    fuzequats.append(float(equat))
    fuzpolar.append(float(polar))         #flips the polar coords
    norm = max(rgb)         #finds max flux value for normalization purposes.
    r, g, b = 255 * redf / norm, 255 * greenf / norm, 255 * bluef / norm
    fuzcolours.append(array([r, g, b]))
    
brightest = max(allfuzaveflux)         #finds the brightest star for normalization purposes

for fuzzy in fuzarray:
    brightness = sqrt(11) * log10((fuzzy[2] / brightest) + 1)            #defines scatter marker 'brightness' based on logarithmic algorithm. Spits out value between 0 and 1
    fuzbright.append(brightness)

fuzcolours = array(fuzcolours)/256        #gets the rgb values between 0 and 1

fig, ax = plt.subplots()
plt.xlabel('Equatorial Angle (deg)')
plt.ylabel('Polar Angle (deg)')
ax.set_facecolor('k')
plt.gca().invert_yaxis()

plt.scatter(fuzequats, fuzpolar, s=fuzbright, c=fuzcolours, marker='.')

figure(figsize=(36,18))             #units are inches
fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
fig.savefig('fuzzymap.png')

totalfuzzy.close()