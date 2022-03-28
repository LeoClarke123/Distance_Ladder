# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 15:07:08 2022

@author: ryanw
"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd
from statistics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def inclination_correction(major, minor, emax, method):
    eccentricity = sqrt(1 - (minor / major)**2)
    if method == "sine":
        velMultiplier = sin(pi/2 * (eccentricity / emax))
    else:
        velMultiplier = cos(pi/2 * (1 - eccentricity / emax))
    return velMultiplier

def radial_distance(galaxX, galaxY, starX, starY):
    deltaX = starX - galaxX
    deltaY = galaxY - starY
    distance = sqrt((deltaX)**2 + (deltaY)**2)
    return distance

dir_path = os.path.dirname(os.path.realpath(__file__))

clusterpath = dir_path.replace("Cluster Analysis", "Data\Camera Images\Star Clusters")

ellipsedata = pd.read_csv(clusterpath.replace("\Star Clusters", "/star clusters.txt"), delimiter=' ')
clusterdistances = pd.read_csv(dir_path.replace("Cluster Analysis", "HR_Diagram/clusterDistances.csv"), delimiter=',')

emax = sqrt(1 - (1.2/7.2)**2)  #maximum eccentricity chosen as edge-on galaxy X245-Y85-N-141

clusterMasses = []
clusterRadii = []
clusterMaxVel = []
clusterGreenLumins = []

StarClusters = []
for clusterFilename in os.listdir(clusterpath):
    clusterName = clusterFilename[0:-4] # remove ".txt" extension - thank you Ciaran!
    StarClusters.append(clusterName)
    
for cluster in StarClusters:
    [clusterX, clusterY, clusterN] = cluster.split("-")
    clusterX = float(clusterX[1:]); clusterY = float(clusterY[1:]); clusterN = float(clusterN[1:])
    
    stardata = pd.read_csv(clusterpath + "/" + cluster + ".txt", delimiter=' ')
    
    starvels = stardata['RadialVelocity'] - median(stardata['RadialVelocity'])  #accounts for the proper motion of the galaxy by subtracting the median velocity from all of the star velocities
    
    major = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,2]; minor = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,3];
    theta = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,4];
    
    
    starvels = abs(starvels * 1 / inclination_correction(major, minor, emax, "sine"))
    
    radius = []
    #the following assumes that *all* galaxies are approximately flat spirals. 
    for index, star in stardata.iterrows():
        radPosition = radial_distance(median(stardata['Equatorial']), median(stardata['Polar']), stardata['Equatorial'].loc[index], stardata['Polar'].loc[index])
        radius.append(radPosition)
        
    distance = clusterdistances.loc[clusterdistances['ClusterName'] == cluster].iloc[0,1]
    distanceMetres = distance * 3.086 * 10**16
    
    clusterRadius = distanceMetres * sin(max(radius) * pi / 180) / 1000
    M = clusterRadius * max(starvels)**2 / (6.637 * 10**-11); clusterMasses.append(M)
    clusterRadii.append(clusterRadius); clusterMaxVel.append(max(starvels))
    GreenLumin = sum(stardata['GreenFlux']) * 4 * pi * distanceMetres**2; clusterGreenLumins.append(GreenLumin)
    
    radius = radius / max(radius)
    
    fig, ax = plt.subplots()
    plt.scatter(radius, starvels, s=10, linewidth=0)
    ax.set_ylabel("Rotational Velocity (km/s)")
    ax.set_xlabel("Radius (prop. of Galactic Radius)")
    ax.set_title(f'Rotation Curve of {cluster}')
    ax.grid(axis='both', which='major')
    fig.savefig(dir_path+f'\\Rotation Curves\\{cluster}.png', dpi=200)
    plt.clf()
    
    


for index, cluster in enumerate(clusterMasses):
    if log10(cluster) > 27:
        del clusterMasses[index]
        del clusterGreenLumins[index]
        del clusterRadii[index]
        del clusterMaxVel[index]

fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterGreenLumins))
ax.set_ylabel("Galaxy V Band Luminosity (log10 W)")
ax.set_xlabel("Galaxy Mass (log10 kg)")
ax.set_title("Galaxy Luminosity vs Mass")

fig.savefig(dir_path + '/Lumin-vs-Mass.png', dpi=200)

    