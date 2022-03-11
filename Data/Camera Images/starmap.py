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
starequats, starpolar, starbright, starcolours, starredshift = [], [], [], [], []
fuzequats, fuzpolar, fuzbright, fuzcolours, fuzredshift = [], [], [], [], []

cm = plt.cm.get_cmap('bwr')     #blue to red colourmap for the redshift graphs. white is the intermediate colour
minvel = 0
maxvel = 0

for row in stardata:
    [name, equat, polar, bluef, greenf, redf, parallax, veloc, distance, period, location] = row.split()
    bluef, greenf, redf, veloc = float(bluef), float(greenf), float(redf), float(veloc)       #cleans up variables
    rgb = [bluef, greenf, redf]
    aveflux = average(rgb)          #averages the r g b flux values
    allstaraveflux.append(aveflux)      #adds it to a list
    stararray.append([equat, polar, aveflux])       
    starequats.append(float(equat))
    starpolar.append(float(polar))         #flips the polar coords
    norm = max(rgb)         #finds max flux value for normalization purposes.
    r, g, b = 255 * redf / norm, 255 * greenf / norm, 255 * bluef / norm            #assigns the rgb values
    starcolours.append(array([r, g, b]))
    #below specifies the new min/max velocity to calibrate the redshift colourbar
    if starredshift != []:
        if veloc < min(starredshift):
            minvel = veloc
        elif veloc > max(starredshift):
            maxvel = veloc
    starredshift.append(veloc)  #adds the velocity of the current star to a redshift array for use in the graph

if maxvel < -minvel:
    maxvel = -minvel
else:
    minvel = -maxvel
    
brightest = max(allstaraveflux)         #finds the brightest star for normalization purposes

for star in stararray:
    brightness = sqrt(11) * log10((star[2] / brightest) + 1)            #defines scatter marker 'brightness' based on logarithmic algorithm. Spits out value between 0 and 1
    starbright.append(brightness)

starcolours = array(starcolours)/256        #gets the rgb values between 0 and 1 (because the scatter function chucked a tantrum)

#below plots the colour star map
fig, ax = plt.subplots()            #initialize axes
ax.set_xlabel('Equatorial Angle (deg)')
ax.set_ylabel('Polar Angle (deg)')
ax.set_facecolor('k')
ax.invert_yaxis()        #flips the y-axis so that polar angle is 0 at due north

plt.scatter(starequats, starpolar, s=starbright, c=starcolours, marker='.')

figure(figsize=(36,18))             #units are inches
fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
fig.savefig(dir_path+'\\starmap.png')

plt.clf()

# make 2 versions, one for analysis (blank name) and one for inclusion in report ('pretty' name)
for (scale,lwidth,name,DPI) in [(0.1,0,'',1200),(5,0,'-pretty',400)]:
    #the below plots the redshift map of the stars
    fig, ax = plt.subplots(figsize=(12,6))  #units are inches
    ax.set_xlabel('Equatorial Angle (deg)')
    ax.set_ylabel('Polar Angle (deg)')
    ax.set_facecolor('k')
    ax.invert_yaxis() 

    # version for analysis
    red = plt.scatter(starequats, starpolar, s=scale, c=starredshift, vmin=minvel, vmax=maxvel, cmap=cm , marker='.', linewidths=lwidth)  #note the colourmap for the redshift amount
    cbar = plt.colorbar(red)
    cbar.set_label('Radial Velocity (km/s)', rotation=90)
    fig.savefig(dir_path+f'\\star-redshift{name}.png', dpi=DPI)


totalpoints.close()
plt.clf()           #clears the current figure in order to create the next one

minvel = 0
maxvel = 0

for row in fuzzydata:       #functionally identical to the star loop
    [name, equat, polar, bluef, greenf, redf, size, veloc, location] = row.split()
    bluef, greenf, redf, veloc = float(bluef), float(greenf), float(redf), float(veloc)
    rgb = [bluef, greenf, redf]
    aveflux = average(rgb)
    allfuzaveflux.append(aveflux)
    fuzarray.append([equat, polar, aveflux])
    fuzequats.append(float(equat))
    fuzpolar.append(float(polar))         #flips the polar coords
    norm = max(rgb)         #finds max flux value for normalization purposes.
    r, g, b = 255 * redf / norm, 255 * greenf / norm, 255 * bluef / norm
    fuzcolours.append(array([r, g, b]))
    fuzredshift.append(veloc)
    if fuzredshift != []:
        if veloc < min(fuzredshift):
            minvel = veloc
    
brightest = max(allfuzaveflux)         #finds the brightest star for normalization purposes

for fuzzy in fuzarray:
    brightness = sqrt(11) * log10((fuzzy[2] / brightest) + 1)            #defines scatter marker 'brightness' based on logarithmic algorithm. Spits out value between 0 and 1
    fuzbright.append(brightness)

fuzcolours = array(fuzcolours)/256        #gets the rgb values between 0 and 1


fig, ax = plt.subplots()
ax.set_xlabel('Equatorial Angle (deg)')
ax.set_ylabel('Polar Angle (deg)')
ax.set_facecolor('k')
ax.invert_yaxis()

plt.scatter(fuzequats, fuzpolar, s=fuzbright, c=fuzcolours, marker='.')

figure(figsize=(36,18))             #units are inches
fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
fig.savefig(dir_path+'\\fuzzymap.png')

plt.clf()

# cm = plt.cm.get_cmap('Blues')

# fig, ax = plt.subplots()
# plt.xlabel('Equatorial Angle (deg)')
# plt.ylabel('Polar Angle (deg)')
# ax.set_facecolor('k')
# plt.gca().invert_yaxis() 

# red = plt.plot(fuzequats, fuzpolar, s=0.1, c=fuzredshift, vmin=-2000, vmax=2000, cmap=cm , marker='.')  #note the colourmap for the redshift amount
# cbar = plt.colorbar(red)
# cbar.set_label('Radial Velocity (km/s)', rotation=90)
# plt.show()

# figure(figsize=(36,18))             #units are inches
# fig.set_dpi(1200)           #sets resolution of image in pixels per square inch?
# fig.savefig(dir_path+'fuzzy-redshift.png')

totalfuzzy.close()