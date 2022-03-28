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
Sclustername = star_c_distances['ClusterName']; Sclusterdist = star_c_distances['Distances']

i = 0
vel = []; Logvel = []
Logphotons = []
GalaxLogDists = []
GalaxDists = []
for source in range(len(xrays)):
    if xrayPhot[source] > 100000:
        for cluster in range(len(starclusters)):
            if in_ellipse(starCx[cluster], starCy[cluster], starChor[cluster], starCvert[cluster], starCangle[cluster], xrayEquat[source], xrayPolar[source]):
                print(xrayPhot[source], "X-Ray Photons in star cluster X =", starCx[cluster], "and Y =",starCy[cluster])
                i += 1
    else:
        for cluster in range(len(galaxclusters)):
            [prefix, xpos, ypos, pop] = GalaxNames[cluster].split('-')
            xpos = float(xpos[-3:]); ypos = float(ypos[-3:]); pop = float(pop[-3:])
            if in_area(xpos, ypos, xrayEquat[source], xrayPolar[source], GalaxVeloc[cluster]):
                #print(xrayPhot[source], " Photons in cluster with speed ", GalaxVeloc[cluster], " with name ", GalaxNames[cluster])
                Logphotons.append(log(xrayPhot[source])); vel.append(GalaxVeloc[cluster]); Logvel.append(-log(abs(GalaxVeloc[cluster])))
                distance = sqrt(max(xrayPhot) * (Sclusterdist.loc[Sclustername == "X187.0-Y122.0-N89"].iloc[0])**2 / xrayPhot[source])
                GalaxLogDists.append(log(distance))
                GalaxDists.append(distance)
                i += 1
                break       #this break prevents double-counting the same X-Ray source for optically close clusters

fig, ax = plt.subplots()            #initialize axes
ax.set_ylabel('Radial Velocity (km/s)')
ax.set_xlabel('Distance from X-Ray Source (pc)')
#ax.invert_xaxis()
plt.scatter(GalaxDists, vel)

z = polyfit(GalaxDists, vel, 1)
p = poly1d(z)

plt.plot(GalaxDists,p(GalaxDists),"r--", linewidth=0.5)
text = f"$y={z[0]:0.5f}\;x{z[1]:+0.3f}$\n$R^2 = {r2_score(vel,p(GalaxDists)):0.3f}$"
plt.gca().text(100000, -800, text)
ax.grid()
ax.ticklabel_format(axis='x', style='sci', scilimits=(0,3), useMathText=True)
ax.set_ylim([-1200,0])
figure(figsize=(20,15))             #units are inches
fig.set_dpi(300)
fig.savefig(dir_path+'\\Velocity vs Distance for Galaxies.png')
