# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 17:01:52 2022

@author: Ryan White     s4499039

Takes a batch of stellar/galactic data and compiles it into neater files, with some analysis included. 
The data format must be in folders in the form:
    DATA
        Back
        .   A
        .   .   01
        .   .   .   fuzzy.txt
        .   .   .   points.txt
        .   .   .   "starname".csv  (optional)
        Up  F   06  ... 

Outputs two files:
    total fuzzy data.txt
    total point-like data.txt
which are exactly as described. 

This program also checks whether stars have associated variable data, and then calculates the periodicity if found data. 

The output data format is in the form of (for fuzzy stars and point-like objects respectively):
    Name    Equatorial Angle (deg)  Polar Angle (deg)   Blue Flux (W/m^2/nm)    Green Flux ("")     Red Flux ("")   Size (arcsec)   Radial Velocity (km/s)  Location (picture)
    Name    Equatorial Angle (deg)  Polar Angle (deg)   Blue Flux (W/m^2/nm)    Green Flux ("")     Red Flux ("")   Parallax (arcsec)   Radial Velocity (km/s)  Distance (pc)   Periodicity (hrs)  Location (picture)

Note that in the decided spherical coordinate system, polar angle is 0 at 'due north' [(0,0) in Up image].
    Positive polar direction is down 
 and equatorial [equat] angle is 0 at the 'front' [(0,0) in the front image, or (0, 90) in polar coords]. 
    Positive equatorial direction is clockwise when viewed from polar north looking down. 
    That is, rotating right is increasing equatorial angle. 

A very big thank you to Dr Benjamin Pope for information /  a tutorial regarding the analysis of variable star data. 
"""

from numpy import *
from astropy.timeseries import LombScargle
from astropy.table import Table
import os 


dir_path = os.path.dirname(os.path.realpath(__file__))  #finds the path of this program to use later

#following initializes all lists needed for program
PNAMES, PEQUATS, PPOLARS, PBLUE, PGREEN, PRED, PPARA, PRVEL, PDIST, PPERIOD, PLOCATION = [], [], [], [], [], [], [], [], [], [], []
FNAMES, FEQUATS, FPOLARS, FBLUE, FGREEN, FRED, FSIZE, FRVEL, FLOCATION = [], [], [], [], [], [], [], [], []

for i in ["Back", "Down", "Front", "Left", "Right", "Up"]:      #each of the six image directions
    for j in ["A", "B", "C", "D", "E", "F"]:                    #each of the image column subdivisions
        for k in ["01", "02", "03", "04", "05", "06"]:          #each of the image row subdivisions
            folder = dir_path + "/" + i + "/" + j + k           #sorts through the data to find data for each subdivision
            fuzzy = open(folder+"/fuzzy.txt", 'r')
            points = open(folder+"/points.txt", 'r')
            variablestars = []          #initialise a list for following loop
            #following for loop checks the folder for the presence of any variable stars, and adds them to a list if found
            for filename in os.listdir(folder):
                if filename.endswith('.csv'):
                    variablestars.append(filename.replace('.csv', ''))
            with fuzzy as f:
                fuzzycontents = f.readlines()[1:]                    #reads the file and puts all content (except line 1) into the list 'contents'
                for row in fuzzycontents:
                    [name, xpos, ypos, bluef, greenf, redf, size, veloc] = row.split()      #splits the row into strings of each data type instead of one long string
                    xpos, ypos, bluef, greenf, redf, size, veloc = float(xpos), float(ypos), float(bluef), float(greenf), float(redf), float(size), float(veloc)
                    #below code converts all data coordinates into a unified, spherical coordinate system
                    if i == "Up":
                        polar = abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Down":
                        polar = 180 - abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Front":
                        polar = abs(ypos - 90)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Back":
                        polar = abs(ypos - 90)
                        equat = 180 + xpos
                    elif i == "Right":
                        polar = abs(ypos - 90)
                        equat = 90 + xpos
                    else:
                        polar = abs(ypos - 90)
                        equat = 270 + xpos
                    #below code appends each variable to the associated string
                    FNAMES.append(name), FEQUATS.append(str(round(equat, 3))), FPOLARS.append(str(round(polar, 3))), FBLUE.append(bluef), FRED.append(redf), FGREEN.append(greenf)
                    FSIZE.append(str(size)), FRVEL.append(veloc), FLOCATION.append(i+j+k)
            fuzzy.close()
            #now to analyse/combine the point data
            with points as p:       #this section's coordinate transform is identical to the fuzzy section, but replaces size variable with parallax angle variable
                pointcontents = p.readlines()[1:]                    #reads the file and puts all content (except line 1) into the list 'pointcontents'
                for row in pointcontents:
                    [name, xpos, ypos, bluef, greenf, redf, parallax, veloc] = row.split()      #splits the row into strings of each data type instead of one long string
                    xpos, ypos, parallax = float(xpos), float(ypos), float(parallax)
                    period = 0          #sets variability period to 0 for stars with no variable nature. Variable stars period is changed in next loop if applicable
                    # following for loop courtesy of Dr Benjamin Pope as per the "Extracting Period Information" section of the PHYS3080 Distance Ladder project info handout
                    for star in variablestars:
                        if name == star:
                            data = Table.read(folder + '/' + name + ".csv", format='ascii')         #reads the .csv data for the variable star
                            time, flux = data['Time'], data['Flux']     #extracts the columns as variables
                            LS = LombScargle(time, flux)            #initializes a Lomb-Scargle algorithm from Astropy
                            freqs = linspace(1/100, 0.45, 10000)        #frequency grid
                            power = LS.power(freqs)     #calculates LS power
                            period = 1 / freqs[argmax(power)]       #finds most likely period from the frequency associated with maximum power
                    
                    #following if statement calculates point distance if parallax angle is not 0. 
                    if parallax > 0.001:
                        distance = round(abs(1 / parallax))
                    else:
                        distance = "undef"
                    
                    #following computes the coordinate transform as before
                    if i == "Up":
                        polar = abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Down":
                        polar = 180 - abs(ypos)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Front":
                        polar = abs(ypos - 90)
                        if xpos < 0:
                            equat = 360 + xpos
                        else:
                            equat = xpos
                    elif i == "Back":
                        polar = abs(ypos - 90)
                        equat = 180 + xpos
                    elif i == "Right":
                        polar = abs(ypos - 90)
                        equat = 90 + xpos
                    else:
                        polar = abs(ypos - 90)
                        equat = 270 + xpos
                    PNAMES.append(name), PEQUATS.append(str(round(equat, 3))), PPOLARS.append(str(round(polar, 3))), PBLUE.append(bluef), PRED.append(redf), PGREEN.append(greenf)
                    PPARA.append(str(parallax)), PRVEL.append(veloc), PDIST.append(distance), PPERIOD.append(str(period)), PLOCATION.append(i+j+k)
            points.close()
print("extra now")

#below code arranges each variable into respective columns for a total data set
totalfuzzy = [FNAMES, FEQUATS, FPOLARS, FBLUE, FGREEN, FRED, FSIZE, FRVEL, FLOCATION]
totalpoints = [PNAMES, PEQUATS, PPOLARS, PBLUE, PGREEN, PRED, PPARA, PRVEL, PDIST, PPERIOD, PLOCATION]

#below code creates names for each column, and writes it to the file
fuzzy = Table(totalfuzzy, names=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Size', 'RadialVelocity', 'Location'])
fuzzy.write('total fuzzy data.txt', overwrite=True, format='ascii')
points = Table(totalpoints, names=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Parallax', 'RadialVelocity', 'Distance', 'Periodicity', 'Location'])  #
points.write('total point-like data.txt', overwrite=True, format='ascii')



