# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 12:15:56 2022

@author: ryanw
"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def in_ellipse(xc, yc, a, b, angle, equat, polar):
    '''
    Finds whether the x/y coords equat/polar are within the ellipse defined by a center (xc, yc), with semi-major axis 'a', semi-minor axis 'b', 
    and angle relative to the positive x-direction 'angle' (in degrees)
    Output is a bool with respect to whether the point is in the ellipse or not. 
    '''
    theta = angle * pi / 180            #converts from degrees to radians
    xdist = equat - xc; ydist = polar - yc      #finds separation from point and ellipse center along each axis
    position = ((xdist * cos(theta) + ydist * sin(theta)) / (a/2))**2 + ((xdist * sin(theta) - ydist * cos(theta)) / (b/2))**2  #finds position of point relative to ellipse
    if position <= 1:       #if within ellipse:
        return True
    else:
        return False
    
def in_area(XC, YC, XS, YS, CVeloc):
    if (XC - 2 <= XS <= XC + 2) and (YC - 2 <= YS <= YC + 2):
        return True
    else:
        return False

dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the file directory

xrays = pd.read_csv(dir_path + "\\Camera Images\\X-Ray Sources.txt", delimiter=' ')
star_c_distances = pd.read_csv(dir_path.replace("\\Data", '') + "\\HR_Diagram\\clusterDistances.csv", delimiter=',')
starclusters = pd.read_csv(dir_path + "\\Camera Images\\star clusters.txt", delimiter=' ')
galaxclusters = pd.read_csv(dir_path + "\\Camera Images\\fuzzy clusters.txt", delimiter=' ')

xrayEquat = xrays['Equatorial']; xrayPolar = xrays['Polar']; xrayPhot = xrays['PhotonsDetected']
starCx = starclusters['x']; starCy = starclusters['y']; starChor = starclusters['horizontal']; starCvert = starclusters['vertical']; starCangle = starclusters['theta']
GalaxNames = galaxclusters['Name']; GalaxNum = galaxclusters['No.Members']; GalaxVeloc = galaxclusters['RadialVelocity']

i = 0
X = []; Xabs = []
Y = []
for source in range(len(xrays)):
    if xrayPhot[source] > 100000:
        for cluster in range(len(starclusters)):
            if in_ellipse(starCx[cluster], starCy[cluster], starChor[cluster], starCvert[cluster], starCangle[cluster], xrayEquat[source], xrayPolar[source]):
                print("h")
                i += 1
                print("i =", i)
    else:
        for cluster in range(len(galaxclusters)):
            [prefix, xpos, ypos, pop] = GalaxNames[cluster].split('-')
            xpos = float(xpos[-3:]); ypos = float(ypos[-3:]); pop = float(pop[-3:])
            if in_area(xpos, ypos, xrayEquat[source], xrayPolar[source], GalaxVeloc[cluster]):
                print("y")
                print(xrayPhot[source], " Photons in cluster with speed ", GalaxVeloc[cluster], " with name ", GalaxNames[cluster])
                Y.append(log(xrayPhot[source])); X.append(GalaxVeloc[cluster]); Xabs.append(log(abs(GalaxVeloc[cluster])))
                i += 1
                print("i =", i)
                break       #this break prevents double-counting the same X-Ray source for optically close clusters

fig, ax = plt.subplots()            #initialize axes
ax.set_xlabel('Log Absolute Radial Velocity (km/s)')
ax.set_ylabel('Log Photons Detected from X-Ray Source')
#ax.invert_xaxis()
plt.scatter(Xabs, Y)
