# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 12:15:56 2022

@author: ryanw

This program takes data in the form of X-Ray Strength/Position, Star Cluster Position/Distance, and Fuzzy Cluster Position to first assign X-Ray sources to fuzzy
clusters and then find the distance to said cluster. 

Distances to distant galaxy clusters are found via the formula:
    
    d_s = sqrt( F_max * (d_max)^2 / F_s )

which takes the distance to the source of maximum flux (d_max and F_max respectively), and the flux of the current x-ray source (F_s) to find the distance to the current source. 
This formula was derived from the inverse square law relation:
    
    intensity1 / intensity2 = distance2^2 / distance1^2

This program outputs a graph that shows the Hubble Constant for this universe, by finding the trendline of a "radial velocity vs distance" relationship. 
"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from sklearn.metrics import r2_score

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
    
def in_area(XC, YC, XS, YS):
    '''Returns a boolean variable as to whether a source seems to be within +/- 2 degrees of a cluster. 
    '''
    if (XC - 2 <= XS <= XC + 2) and (YC - 2 <= YS <= YC + 2):
        return True
    else:
        return False

dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the file directory

#read the data from corresponding files
xrays = pd.read_csv(dir_path + "\\Camera Images\\X-Ray Sources.txt", delimiter=' ')
star_c_distances = pd.read_csv(dir_path.replace("\\Data", '') + "\\HR_Diagram\\clusterDistances.csv", delimiter=',')
starclusters = pd.read_csv(dir_path + "\\Camera Images\\star clusters.txt", delimiter=' ')
galaxclusters = pd.read_csv(dir_path + "\\Camera Images\\fuzzy clusters.txt", delimiter=' ')

#give specific columns some easier-to-call names
xrayEquat = xrays['Equatorial']; xrayPolar = xrays['Polar']; xrayPhot = xrays['PhotonsDetected']
starCx = starclusters['x']; starCy = starclusters['y']; starChor = starclusters['horizontal']; starCvert = starclusters['vertical']; starCangle = starclusters['theta']
GalaxNames = galaxclusters['Name']; GalaxNum = galaxclusters['No.Members']; GalaxVeloc = galaxclusters['RadialVelocity']
Sclustername = star_c_distances['ClusterName']; Sclusterdist = star_c_distances['Distances']

i = 0

#initialise required lists to work with
vel = []; Logvel = []
Logphotons = []
GalaxLogDists = []
GalaxDists = []

#now to iterate over each X-Ray source!
for source in range(len(xrays)):
    if xrayPhot[source] > 100000:       #this must be a very close X-Ray source and would be within a galaxy with resolved stars
        for cluster in range(len(starclusters)):
            if in_ellipse(starCx[cluster], starCy[cluster], starChor[cluster], starCvert[cluster], starCangle[cluster], xrayEquat[source], xrayPolar[source]):
                print(xrayPhot[source], "X-Ray Photons in star cluster X =", starCx[cluster], "and Y =",starCy[cluster])
                i += 1
    else:       #this is for more distant, unresolved galaxies
        for cluster in range(len(galaxclusters)):       #for loop iterating over all of the galaxy clusters
            [prefix, xpos, ypos, pop] = GalaxNames[cluster].split('-')          #gets data from cluster name
            xpos = float(xpos[-3:]); ypos = float(ypos[-3:]); pop = float(pop[-3:])     #cleans up data from the name of the cluster
            if in_area(xpos, ypos, xrayEquat[source], xrayPolar[source]):       #checks if cluster aligns optically with a galaxy cluster
                #print(xrayPhot[source], " Photons in cluster with speed ", GalaxVeloc[cluster], " with name ", GalaxNames[cluster])
                
                distance = sqrt(max(xrayPhot) * (Sclusterdist.loc[Sclustername == "X187.0-Y122.0-N89"].iloc[0])**2 / xrayPhot[source]) #calculates distance to this cluster
                #now to append the data to the corresponding lists
                Logphotons.append(log(xrayPhot[source])); vel.append(GalaxVeloc[cluster]); Logvel.append(-log(abs(GalaxVeloc[cluster])))
                GalaxLogDists.append(log(distance)); GalaxDists.append(distance)
                i += 1
                break       #this break prevents double-counting the same X-Ray source for optically close clusters

#the following code creates a velocity-distance graph for the distant galaxy clusters. 
fig, ax = plt.subplots()            #initialize axes
ax.set_ylabel('Radial Velocity (km/s)')
ax.set_xlabel('Distance from X-Ray Source (pc)')
#ax.invert_xaxis()
plt.scatter(GalaxDists, vel)

z = polyfit(GalaxDists, vel, 1)     #this finds the linear fit for the data
p = poly1d(z)

plt.plot(GalaxDists,p(GalaxDists),"r--", linewidth=0.5)     #plot the trendline on top of the data

text = f"$y={z[0]:0.5f}\;x{z[1]:+0.3f}$\n$R^2 = {r2_score(vel,p(GalaxDists)):0.3f}$"        #this defines the text to add onto the graph
plt.gca().text(100000, -800, text)      #chooses the location for the text on the graph
ax.grid()
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,3), useMathText=True)       #scientific notation for the x-axis
ax.set_ylim([-1200,0])
figure(figsize=(20,15))             #units are inches
fig.set_dpi(300)
fig.savefig(dir_path+'\\Velocity vs Distance for Galaxies.png')
