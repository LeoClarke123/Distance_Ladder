# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 17:01:52 2022

@author: Ryan White     s4499039

Takes a batch of stellar/galactic/xray data and compiles it into neater files, with some analysis included and a transform of coordinates from the 6 images to a unified system.
The data format must be in folders in the form:
    DATA
        Back
        .   A
        .   .   01
        .   .   .   fuzzy.txt
        .   .   .   points.txt
        .   .   .   "starname".csv  (optional)
        Up  F   06  ... 
with the exception of the X-Ray data file, which must be in the directory *above* the one this file is in. 

Outputs four files:
    total fuzzy data.txt
    total point-like data.txt
    calibrated star data.txt
    X-Ray Sources.txt
which are exactly as described, where "calibrated star data.txt" is a list of stars for which the distance was inferred by parallax (important for further analysis). 

This program also checks whether stars have associated variable data, and then calculates the periodicity if data is found. 

The output data format is in the form of (for fuzzy stars and point-like/calibrated objects respectively):
    Name    Equatorial Angle (deg)  Polar Angle (deg)   Blue Flux (W/m^2/nm)    Green Flux ("")     Red Flux ("")   Size (arcsec)   Radial Velocity (km/s)  Distance (pc)   Location (picture)
    Name    Equatorial Angle (deg)  Polar Angle (deg)   Blue Flux (W/m^2/nm)    Green Flux ("")     Red Flux ("")   Parallax (arcsec)   Radial Velocity (km/s)  Distance (pc)   Periodicity (hrs)  Period Uncertainty (+/- hrs)  Location (picture)

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
import pandas as pd 
from errorFormulae import period_uncert

def coord_transform(direction, xpos, ypos):
    '''
    It's gotten to the point that even I'm struggling to follow what this does. For the Up and Down directions, a square to circle coordinate transform is done.
    
    '''
    if direction == "Up":
        u = abs(xpos) * sqrt(1 - (1/2 * (ypos/45)**2))
        v = abs(ypos) * sqrt(1 - (1/2 * (xpos/45)**2))
        polar = sqrt(u**2 + v**2)
        if xpos > 0 and ypos < 0:
            equat = degrees(arctan(xpos / abs(ypos)))  
        elif xpos > 0 and ypos > 0:
            equat = 180 - degrees(arctan(xpos / ypos))
        elif xpos < 0 and ypos > 0:
            equat = 180 + degrees(arctan(abs(xpos) / ypos))
        else:
            equat = 360 - degrees(arctan(abs(xpos) / abs(ypos)))
    elif direction == "Down":
        u = abs(xpos) * sqrt(1 - (1/2 * (ypos/45)**2))
        v = abs(ypos) * sqrt(1 - (1/2 * (xpos/45)**2))
        polar = 180 - sqrt(u**2 + v**2)
        if xpos > 0 and ypos > 0:
            equat = degrees(arctan(xpos / ypos))   
        elif xpos > 0 and ypos < 0:
            equat = 180 - degrees(arctan(xpos / abs(ypos)))
        elif xpos < 0 and ypos < 0:
            equat = 180 + degrees(arctan(abs(xpos) / abs(ypos)))
        else:
            equat = 360 - degrees(arctan(abs(xpos) / ypos))
    elif direction == "Front":
        polar = abs(ypos - 90)
        if xpos < 0:
            equat = 360 + xpos
        else:
            equat = xpos
    elif direction == "Back":
        polar = abs(ypos - 90)
        equat = 180 + xpos
    elif direction == "Right":
        polar = abs(ypos - 90)
        equat = 90 + xpos
    else:
        polar = abs(ypos - 90)
        equat = 270 + xpos
    return [equat, polar]


dir_path = os.path.dirname(os.path.realpath(__file__))  #finds the path of this program to use later

#following initializes all lists needed for program; P for point-like objects, C for calibrated stars, and F for fuzzy objects
PNAMES, PEQUATS, PPOLARS, PBLUE, PGREEN, PRED, PPARA, PRVEL, PDIST, PPERIOD, PPERIODUNC, PLOCATION = [], [], [], [], [], [], [], [], [], [], [], []
CNAMES, CEQUATS, CPOLARS, CBLUE, CGREEN, CRED, CPARA, CRVEL, CDIST, CPERIOD, CPERIODUNC, CLOCATION = [], [], [], [], [], [], [], [], [], [], [], []
FNAMES, FEQUATS, FPOLARS, FBLUE, FGREEN, FRED, FSIZE, FRVEL, FDIST, FLOCATION = [], [], [], [], [], [], [], [], [], []

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
                    [equat, polar] = coord_transform(i, xpos, ypos)
                    #below code appends each variable to the associated string
                    
                    distance = (veloc - 1.643) / -0.00263       #this function derived from the x-ray flash to galaxy cluster distance/velocity relation
                    
                    FNAMES.append(name), FEQUATS.append(str(round(equat, 3))), FPOLARS.append(str(round(polar, 3))), FBLUE.append(bluef), FRED.append(redf), FGREEN.append(greenf)
                    FSIZE.append(str(size)), FRVEL.append(veloc), FLOCATION.append(i+j+k), FDIST.append(distance)
            fuzzy.close()
            #now to analyse/combine the point data
            with points as p:       #this section's coordinate transform is identical to the fuzzy section, but replaces size variable with parallax angle variable
                pointcontents = p.readlines()[1:]                    #reads the file and puts all content (except line 1) into the list 'pointcontents'
                for row in pointcontents:
                    [name, xpos, ypos, bluef, greenf, redf, parallax, veloc] = row.split()      #splits the row into strings of each data type instead of one long string
                    xpos, ypos, parallax = float(xpos), float(ypos), float(parallax)
                    period = 0          #sets variability period to 0 for stars with no variable nature. Variable stars period is changed in next loop if applicable
                    period_unc = 0
                    # following for loop courtesy of Dr Benjamin Pope as per the "Extracting Period Information" section of the PHYS3080 Distance Ladder project info handout
                    for star in variablestars:
                        if name == star:
                            data = Table.read(folder + '/' + name + ".csv", format='ascii')         #reads the .csv data for the variable star
                            time, flux = data['Time'], data['Flux']     #extracts the columns as variables
                            LS = LombScargle(time, flux)            #initializes a Lomb-Scargle algorithm from Astropy
                            freqs = linspace(1/100, 0.45, 10000)        #frequency grid
                            power = LS.power(freqs)     #calculates LS power
                            period = 1 / freqs[argmax(power)]       #finds most likely period from the frequency associated with maximum power
                            
                            period_unc = period_uncert(time, flux, runs = 20,  yerr = 0.015)    #finds the uncertainty in period from monte carlo sim developed by Leo Clarke
                            print("done")
                    #following computes the coordinate transform as before
                    [equat, polar] = coord_transform(i, xpos, ypos)
                        
                    #following if statement calculates point distance if parallax angle is greater than some threshold. This also appends the calibrated lists with the star data
                    if parallax > 0.005:
                        distance = round(abs(1 / parallax))
                        CNAMES.append(name), CEQUATS.append(str(round(equat, 3))), CPOLARS.append(str(round(polar, 3))), CBLUE.append(bluef), CRED.append(redf), CGREEN.append(greenf)
                        CPARA.append(str(parallax)), CRVEL.append(veloc), CDIST.append(distance), CPERIOD.append(str(period)), CPERIODUNC.append(str(period_unc)), CLOCATION.append(i+j+k)
                    else:
                        distance = "undef"
                    PNAMES.append(name), PEQUATS.append(str(round(equat, 3))), PPOLARS.append(str(round(polar, 3))), PBLUE.append(bluef), PRED.append(redf), PGREEN.append(greenf)
                    PPARA.append(str(parallax)), PRVEL.append(veloc), PDIST.append(distance), PPERIOD.append(str(period)), PPERIODUNC.append(str(period_unc)), PLOCATION.append(i+j+k)
            points.close()

#below code arranges each variable into respective columns for a total data set
totalfuzzy = [FNAMES, FEQUATS, FPOLARS, FBLUE, FGREEN, FRED, FSIZE, FRVEL, FDIST, FLOCATION]
totalpoints = [PNAMES, PEQUATS, PPOLARS, PBLUE, PGREEN, PRED, PPARA, PRVEL, PDIST, PPERIOD, PPERIODUNC, PLOCATION]
totalcalibrated = [CNAMES, CEQUATS, CPOLARS, CBLUE, CGREEN, CRED, CPARA, CRVEL, CDIST, CPERIOD, CPERIODUNC, CLOCATION]

#below code creates names for each column, and writes it to the file
fuzzy = Table(totalfuzzy, names=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Size', 'RadialVelocity', 'Distance', 'Location'])
fuzzy.write('total fuzzy data.txt', overwrite=True, format='ascii')

points = Table(totalpoints, names=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Parallax', 'RadialVelocity', 'Distance', 'Periodicity', 'PeriodUnc', 'Location'])  #
points.write('total point-like data.txt', overwrite=True, format='ascii')

calibrated = Table(totalcalibrated, names=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Parallax', 'RadialVelocity', 'Distance', 'Periodicity', 'PeriodUnc', 'Location'])
calibrated.write('calibrated star data.txt', overwrite=True, format='ascii')

#now to transform coordinates of the X-Ray Data
xrays = pd.read_csv(dir_path.replace("Camera Images", '') + 'X-Ray All Sky Data.txt', delimiter=' ')
xrayDirection = xrays['CameraDirection']; xrayX = xrays['X']; xrayY = xrays['Y']; xrayPhotons = xrays['PhotonsDetected']

xraydata = []
for source in range(len(xrays)):
    xpos = xrayX[source]; ypos = xrayY[source] #initialize variables
    
    [equat, polar] = coord_transform(xrayDirection[source], xpos, ypos) #transform coordinates
    
    xraydata.append([equat, polar, xrayPhotons[source]]) #add new coord data to list

NewXRayFile = pd.DataFrame(xraydata, columns=['Equatorial', 'Polar', 'PhotonsDetected'])
NewXRayFile.to_csv(dir_path+'/X-Ray Sources.txt', index=None, sep=' ')    #writes xray data to file