# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 15:07:08 2022

@author: ryanw

This program analyses nearby (and resolved) galaxies to produce a rough rotation curve for each. It also does some analysis on these galaxies to infer key relationships between
galaxy properties (to mixed results). 

This program currently outputs a rotation curve for each nearby galaxy, as well as Luminosity-Mass graph. 
"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd
from statistics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def inclination_correction(major, minor, emax, method):
    '''This function increases the rotation velocity of observed stars in galaxies depending on some inferred inclination. 
    It does this by analysing the eccentricity of the observed galaxy shape, and comparing it to a near-perfect 'edge-on' galaxy, with eccentricity given by emax.
    It will calculate the eccentricity of the observed galaxy by it's major and minor axes, and then use sine/cosine (depending on what the user wants), to infer
    a multiplier for which, when the observed rotational velocity is multiplied by it, will return an estimation of the *actual* rotational velocity. 
    '''
    eccentricity = sqrt(1 - (minor / major)**2)
    if method == "sine":
        velMultiplier = sin(pi/2 * (eccentricity / emax))
    else:
        velMultiplier = cos(pi/2 * (1 - eccentricity / emax))
    return velMultiplier

def radial_distance(galaxX, galaxY, starX, starY):
    '''This function finds an approximate distance that a star is from the center of the galaxy.
    (Does not account for inclination)
    '''
    deltaX = starX - galaxX
    deltaY = galaxY - starY
    distance = sqrt((deltaX)**2 + (deltaY)**2)
    return distance

dir_path = os.path.dirname(os.path.realpath(__file__))

clusterpath = dir_path.replace("Cluster Analysis", "Data\Camera Images\Star Clusters")
galaxData = pd.read_csv(clusterpath.replace("\Star Clusters", "/total fuzzy data.txt"), delimiter=' ')
ellipsedata = pd.read_csv(clusterpath.replace("\Star Clusters", "/star clusters.txt"), delimiter=' ')
clusterdistances = pd.read_csv(dir_path.replace("Cluster Analysis", "HR_Diagram/clusterDistances.csv"), delimiter=',')

emax = sqrt(1 - (1.2/7.2)**2)  #maximum eccentricity chosen as edge-on galaxy X245-Y85-N-141

clusterMasses = []
clusterRadii = []
clusterMaxVel = []
clusterGreenLumins = []

StarClusters = []
#the following loop adds all clusternames to a list StarClusters
for clusterFilename in os.listdir(clusterpath):
    clusterName = clusterFilename[0:-4] # remove ".txt" extension - thank you Ciaran!
    StarClusters.append(clusterName)
    
#now to analyse each cluster
for cluster in StarClusters:
    [clusterX, clusterY, clusterN] = cluster.split("-")
    clusterX = float(clusterX[1:]); clusterY = float(clusterY[1:]); clusterN = float(clusterN[1:])      #cleans up the data coming from the cluster name
    
    stardata = pd.read_csv(clusterpath + "/" + cluster + ".txt", delimiter=' ')     #loads the cluster data
    
    starvels = stardata['RadialVelocity'] - median(stardata['RadialVelocity'])  #accounts for the proper motion of the galaxy by subtracting the median velocity from all of the star velocities
    
    #finds the observational characteristics of the nearby galaxy
    major = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,2]
    minor = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,3]
    theta = ellipsedata.loc[ellipsedata['x'] == clusterX].iloc[0,4]
    
    
    starvels = abs(starvels * 1 / inclination_correction(major, minor, emax, "cos"))   #estimates the actual rotation velocity via an inclination correction
    
    radius = []
    #the following loop estimates the radial position that a star is from the center of the galaxy
    #the following assumes that *all* galaxies are approximately flat spirals. 
    for index, star in stardata.iterrows():
        radPosition = radial_distance(median(stardata['Equatorial']), median(stardata['Polar']), stardata['Equatorial'].loc[index], stardata['Polar'].loc[index])
        radius.append(radPosition)
    
    #the following finds the distance to the nearby galaxy
    distance = clusterdistances.loc[clusterdistances['ClusterName'] == cluster].iloc[0,1]
    distanceMetres = distance * 3.086 * 10**16
    
    
    clusterRadius = distanceMetres * sin(max(radius) * pi / 180) / 1000     #estimates galaxy radius from it's optical extent and the distance to it
    M = clusterRadius * max(starvels)**2 / (6.637 * 10**-11); clusterMasses.append(M)       #estimates TOTAL galaxy mass from the virial theorem
    clusterRadii.append(clusterRadius); clusterMaxVel.append(max(starvels))     #add data to lists
    #the following finds galaxy luminosity by adding the luminosity of all of the stars
    GreenLumin = sum(stardata['GreenFlux']) * 4 * pi * distanceMetres**2; 
    clusterGreenLumins.append(GreenLumin)     #append luminosity to list of all galaxy luminosities
    
    #the following plots a rotation curve of the currently analysed galaxy
    radius = radius / max(radius)       #finds proportion of galactic radius that star is away from the center
    fig, ax = plt.subplots()
    plt.scatter(radius, starvels, s=10, linewidth=0)
    ax.set_ylabel("Rotational Velocity (km/s)")
    ax.set_xlabel("Radius (prop. of Galactic Radius)")
    #ax.set_title(f'Rotation Curve of {cluster}')
    ax.grid(axis='both', which='major')
    fig.savefig(dir_path+f'\\Rotation Curves\\{cluster}.png', dpi=200)
    plt.clf()
    
    

#the following removes severe outliers from the data
for index, cluster in enumerate(clusterMasses):
    if log10(cluster) > 27:
        del clusterMasses[index]
        del clusterGreenLumins[index]
        del clusterRadii[index]
        del clusterMaxVel[index]

#this plots luminosity vs galaxy mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterGreenLumins))
ax.set_ylabel("Log10 Galaxy V Band Luminosity (W)")
ax.set_xlabel("Log10 Galaxy Mass (kg)")
#ax.set_title("Galaxy Luminosity vs Mass")

fig.savefig(dir_path + '/Lumin-vs-Mass.png', dpi=200)

plt.clf()

#this calculates and plots the colour of a galaxy based on distance
fig, ax = plt.subplots()
colour = log(galaxData['RedFlux'] / galaxData['BlueFlux'])
plt.scatter(log10(galaxData['Distance']), colour, s=0.1)
ax.set_ylabel("Galaxy Colour ($M_B - M_R$)")
ax.set_xlabel("Galaxy Distance (pc)")
#ax.set_title("Galaxy Colour vs Distance")

fig.savefig(dir_path + '/Colour-vs-Distance.png', dpi=200)

plt.clf()

#this plots the galaxy radii vs the mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterRadii))
ax.set_ylabel("Log10 Galaxy Radius (metres)")
ax.set_xlabel("Log10 Galaxy Mass (kg)")
#ax.set_title("Galaxy Radius vs Mass")

fig.savefig(dir_path + '/Radius vs Mass.png', dpi=200)
    
plt.clf()

#this plots the cluster rotational velocity vs mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterMaxVel))
ax.set_ylabel("Log10 Galaxy Rot Velocity (km/s)")
ax.set_xlabel("Log10 Galaxy Mass (kg)")
#ax.set_title("Galaxy Rot Velocity vs Mass")

fig.savefig(dir_path + '/Rot Velocity vs Mass.png', dpi=200)

plt.clf()

#this plots the size of a galaxy vs its distance
fig, ax = plt.subplots()
plt.scatter(log10((galaxData['RadialVelocity'] - 1.643) / -0.00263), log10(((galaxData['RadialVelocity'] - 1.643) / -0.00263) * tan(galaxData['Size'])), s=0.1)
ax.set_ylabel("Log10 Galaxy Size (pc)")
ax.set_xlabel("Log10 Galaxy Distance (pc)")
#ax.set_title("Galaxy Size vs Distance")

fig.savefig(dir_path + '/Size vs Distance.png', dpi=200)