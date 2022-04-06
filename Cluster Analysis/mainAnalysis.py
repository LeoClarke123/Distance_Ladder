# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 15:07:08 2022

@author: Ryan White     s4499039

This program analyses nearby (and resolved) galaxies to produce a rough rotation curve for each. It also does some analysis on these galaxies to infer key relationships between
galaxy properties (to mixed results). 

This program currently outputs a rotation curve for each nearby galaxy, as well as a host of graphs about galaxy properties. 

The algorithm to produce approximately equal area cells for a rectangle mapping of a sphere was adapted from:
Malkin, Z, 2019. "A NEW EQUAL-AREA ISOLATITUDINAL GRID ON A SPHERICAL SURFACE", https://arxiv.org/pdf/1909.04701.pdf


"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd
from statistics import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.patches as patches
import matplotlib.colors as colours
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from array import *

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
galaxClusters = pd.read_csv(clusterpath.replace("\Star Clusters", "/fuzzy clusters.txt"), delimiter=' ')

emax = sqrt(1 - (1.2/7.2)**2)  #maximum eccentricity chosen as edge-on galaxy X245-Y85-N-141

clusterMasses = []
clusterRadii = []
clusterMaxVel = []
clusterGreenLumins = []
radUnc = []; massUnc = []; logMassUnc = []
allTrends = []

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
    
    radius = []; radiusUnc = [];
    #the following loop estimates the radial position that a star is from the center of the galaxy
    #the following assumes that *all* galaxies are approximately flat spirals. 
    for index, star in stardata.iterrows():
        radPosition = radial_distance(median(stardata['Equatorial']), median(stardata['Polar']), stardata['Equatorial'].loc[index], stardata['Polar'].loc[index])
        radius.append(radPosition)
        radiusUnc.append(1/2 * radPosition * (2*sin(0.2 * pi/180)**2))
    
    
    #the following finds the distance to the nearby galaxy
    distance = clusterdistances.loc[clusterdistances['ClusterName'] == cluster].iloc[0,1]
    distanceMetres = distance * 3.086 * 10**16
    
    clusterRadius = distanceMetres * sin(max(radius) * pi / 180) / 1000     #estimates galaxy radius from it's optical extent and the distance to it
    radiusUnc = 0.434 * (sin(0.2 * pi/180) * distance / (clusterRadius / (3.086 * 10**13)))
    radUnc.append(radiusUnc)
    M = clusterRadius * (max(starvels) * 1000)**2 / (6.67 * 10**-11); clusterMasses.append(M)       #estimates TOTAL galaxy mass from the virial theorem
    clusterRadii.append(clusterRadius / (3.086 * 10**13)); clusterMaxVel.append(max(starvels))     #add data to lists
    
    #the following loop implements a monte carlo estimation of the galaxy's mass uncertainty
    i=0
    runs = 20
    massUncDist = []
    while i < runs:      
        radius_jitter = clusterRadius + (radiusUnc * random.uniform(low=-1, high=1, size=len(StarClusters)))
        vel_jitter = (max(starvels) * 1000) + ((max(starvels) * 1000) * 0.03 * random.uniform(low=-1, high=1, size=len(StarClusters)))
        massUncDist.append(radius_jitter * vel_jitter**2 / (6.67 * 10**-11))
        i += 1
    massUnc.append(std(massUncDist))
    logMassUnc.append(0.464 * std(massUncDist) / M)
    
    #the following finds galaxy luminosity by adding the luminosity of all of the stars
    GreenLumin = sum(stardata['GreenFlux']) * 4 * pi * distanceMetres**2; 
    clusterGreenLumins.append(GreenLumin)     #append luminosity to list of all galaxy luminosities
    
    #the following plots a rotation curve of the currently analysed galaxy
    radius = radius / max(radius)       #finds proportion of galactic radius that star is away from the center
    fig, ax = plt.subplots()
    plt.scatter(radius, starvels, s=5, linewidth=0)
    plt.errorbar(radius, starvels, xerr=radiusUnc, yerr=0.03, fmt=',', linewidth=0.7)
    ax.set_ylabel("Rotational Velocity (km s$^{-1}$)")
    ax.set_xlabel("Radius (prop. of Galactic Radius)")
    #ax.set_title(f'Rotation Curve of {cluster}')
    ax.grid(axis='both', which='major')
    
    #the below few loops finds the "top" of the rotation curve, since that most accurately represents the *actual* rotation curve
    indRadVel = []
    for index, velocity in enumerate(starvels):
        indRadVel.append([index, radius[index], velocity])
    SortedRadVels = sorted(indRadVel, key=lambda x: x[1])
    trendRads = []; trendVels = []
    for index, value in enumerate(SortedRadVels):
        ind, rad, vel = value
        if index == 0:
            trendRads.append(rad)
            trendVels.append(vel)
        elif vel > max(trendVels):
            trendRads.append(rad)
            trendVels.append(vel)
    
    #the below attempts to plot a trendline of the rotation curve
    x = arange(min(radius), max(radius), 0.01)      #x bounds for the trendline
    z,cov = polyfit(log10(trendRads), trendVels, 1, cov=True)     #this finds the logarithmic fit for the data
    p = poly1d(z)
    gradUnc, intUnc = sqrt(diag(cov))
    upper = (z[0] + gradUnc) * log10(x) + (z[1] - intUnc)
    lower = (z[0] - gradUnc) * log10(x) + (z[1] + intUnc)
    plt.plot(x, p(log10(x)),"r--", linewidth=0.5)
    plt.fill_between(x, lower, upper, color='r', alpha=0.2)
    plt.ylim([0, (1.1) * max(starvels)])
    
    allTrends.append(p / p(max(x)))
    
    fig.savefig(dir_path+f'\\Rotation Curves\\{cluster}.png', dpi=200)
    plt.close(fig)
    plt.clf()
    
x = arange(0.01, 1.01, 0.01)
fig, ax = plt.subplots()
for trendline in allTrends:
    plt.plot(x, trendline(log10(x)), linewidth=0.5)
ax.set_ylabel("Rotational Velocity (prop. of v$_{max}$)")
ax.set_xlabel("Radius (prop. of Galactic Radius)")
plt.ylim([0, 0.8])
fig.savefig(dir_path+'\\Rotation Curves\\ALLCurves.png', dpi=300)

plt.clf()
    

#the following removes severe outliers from the data
for index, cluster in enumerate(clusterMasses):
    if log10(cluster) > 33.5:
        del clusterMasses[index]
        del clusterGreenLumins[index]
        del clusterRadii[index]
        del clusterMaxVel[index]
        del radUnc[index]
        del massUnc[index]
        del logMassUnc[index]


clusterCoords = []
clusterNames = galaxClusters['Name']
for galaxyCluster in clusterNames:
    prefix, xtext, ytext, ntext = galaxyCluster.split('-')
    xtext = float(xtext[-3:]); ytext = float(ytext[-3:])
    clusterCoords.append((xtext, ytext))

rings = 14      #number of rings to divide rectangle into. make sure that this number is even
height = 180/rings      #height of each ring (in degrees)
ringArea = 4 * pi / rings       #each ring should have same area 
upper = pi / 2      #starts with the upper bound of the graph at pi/2 radians (north pole)
ic = 0          #internal counter to check how many galaxy clusters the program identifies overall
scaling = 0.05      #sets size of graphs as a proprtion of 360x180 inches

fig, ax = plt.subplots(figsize=(scaling*360, scaling*180))
cmap = plt.cm.get_cmap('summer')
rects = []; antirects = []; totCellCounts = []  #initialise lists to store data in, to then graph later

for ring in range(1, int((rings / 2) + 1)):       #iterates over the northern hemisphere's rings, with an "antipolar" calculation done for the southern hemisphere
    centralLat = upper - ((height / 2) * pi / 180)      #finds the central latitude (radians) of the current ring
    dLi = height * (1 / cos(centralLat))        #chooses how many horizontal subdivisions of the ring to apply, with more divisions closer to the equator
    lower = upper - (height * pi/180)           #since all rings are equal in width, this finds the lower bound of the ring just based on the ring height and upper bound
    PolarUpper = abs((180 / pi) * (upper - pi/2))   #converts the upper bound in radians to the desired coordinate system (degrees, with 0 deg at north pole)
    PolarLower = PolarUpper + height    #as above, but for the lower bound
    AntiPolarU = 180 - PolarLower       #finds upper bound of southern hemisphere cell
    AntiPolarL = AntiPolarU + height    #as above, but lower bound
    cells = int(round(360/dLi))     #finds an integer number of cells to subdivide the equatorial axis into
    
    for i in range(0, cells):       #now to check how many galaxy clusters are in each cell in this ring
        cellLeft = i * 360 / cells      #finds left coord of the current cell
        cellRight = cellLeft + 360 / cells          #finds right coord of the current cell
        CellCount = 0; AntiCellCount = 0;       #initializes variables for the number of clusters in each cell
        for cluster in clusterCoords:       #iterates over all galaxy clusters to check if it's in this cell
            x, y = cluster
            if (cellLeft <= x <= cellRight) and (PolarUpper <= y <= PolarLower):        #if cluster is in this cell, add it to northern hemisphere counter
                CellCount += 1
                ic += 1
            elif (cellLeft <= x <= cellRight) and (AntiPolarU <= y <= AntiPolarL):  #if cluster is in the anti-cell, add it to the southern hemisphere counter
                AntiCellCount += 1
                ic += 1
        totCellCounts.append(CellCount); totCellCounts.append(AntiCellCount)        #add counts to a list (for colour normalization later on)
        Rwidth = cellRight - cellLeft       #width of rectangle cell
        Rheight = PolarUpper - PolarLower   #height of rectangle cell
        rects.append((cellLeft, PolarLower, Rwidth, Rheight, CellCount))        #add the rectangle properties to a list to plot later
        antirects.append((cellLeft, AntiPolarL, Rwidth, Rheight, AntiCellCount))        #as above, but for southern hemisphere rectangles   
    upper = lower       #changes upper bound for next ring into lower bound of current ring

ColourNorm = colours.Normalize(vmin = 0, vmax = max(totCellCounts))     #normalize the colour map across the cell populations

for rect in rects:      #for all cells in northern hemisphere
    cellLeft, PolarLower, Rwidth, Rheight, CellCount = rect             #gets rectangle properties from the list
    rectangle = patches.Rectangle(xy=(cellLeft, PolarLower), width=Rwidth, height=Rheight, edgecolor='k', facecolor=cmap(ColourNorm(CellCount)), fill=True)
    ax.add_patch(rectangle)         #plots the rectangle on the figure
for rect in antirects:      #as above but for the southern hemisphere
    cellLeft, PolarLower, Rwidth, ARheight, CellCount = rect
    rectangle = patches.Rectangle(xy=(cellLeft, PolarLower), width=Rwidth, height=ARheight, edgecolor='k', facecolor=cmap(ColourNorm(CellCount)), fill=True)
    ax.add_patch(rectangle)

ax.set_xlabel("Equatorial Angle (deg)")
ax.set_ylabel("Polar Angle (deg)")
plt.xlim([0,360])
plt.ylim([0,180])
ax.invert_yaxis()       #flip the axis to match the coordinate system

sm = plt.cm.ScalarMappable(cmap=cmap, norm=ColourNorm)
sm.set_array([])
cbar = plt.colorbar(sm)
cbar.set_label('Cluster Population in Cell', rotation=90)

fig.savefig(dir_path + '/Isotropy.png', dpi=200)    
plt.clf()

print("The mean cluster population per cell is", mean(totCellCounts), "with a standard deviation of", std(totCellCounts))
fig, ax = plt.subplots()
plt.hist(totCellCounts, bins=arange(0,max(totCellCounts)+2,1)-0.5, edgecolor='k')
ax.set_xlabel("Galaxy Clusters in Cell")
ax.set_ylabel("Number of Cells")

fig.savefig(dir_path + '/Cluster Cell Count Histogram.png', dpi=200)


#this plots luminosity vs galaxy mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterGreenLumins), s=3)
ax.set_ylabel("log$_{10}$ Galaxy V Band Luminosity (W m$^{-2}$ nm$^{-1}$)")
ax.set_xlabel("log$_{10}$ Galaxy Mass (kg)")
# x = arange(min(log10(clusterMasses)), max(log10(clusterMasses)), 0.1)
# z,cov = polyfit(log10(clusterMasses), log10(clusterGreenLumins), 1, cov=True)     #this finds the linear fit for the data
# p = poly1d(z)
# gradUnc, intUnc = sqrt(diag(cov))
# upper = (z[0] + gradUnc) * x + (z[1] - intUnc)
# lower = (z[0] - gradUnc) * x + (z[1] + intUnc)
# plt.plot(log10(clusterMasses),p(log10(clusterMasses)),"r--", linewidth=0.5)
# plt.fill_between(x, lower, upper, color='r', alpha=0.2)

# print(f"log10(L_G)=({z[0]:0.5f} \pm {round(gradUnc,5)}) M + ({z[1]:0.3f} \pm {round(intUnc, 2)}) km/s \nR^2 = {r2_score(log10(clusterGreenLumins),p(log10(clusterMasses))):0.3f}")
plt.errorbar(log10(clusterMasses), log10(clusterGreenLumins), xerr=logMassUnc, yerr=0.03, fmt=',', linewidth=0.5)
#ax.set_title("Galaxy Luminosity vs Mass")

fig.savefig(dir_path + '/Lumin-vs-Mass.png', dpi=200)

plt.clf()

#this calculates and plots the colour of a galaxy based on distance
fig, ax = plt.subplots()
colour = log(galaxData['RedFlux'] / galaxData['BlueFlux'])
plt.scatter(log10(galaxData['Distance']), colour, s=0.1)
ax.set_ylabel("Galaxy Colour ($M_B - M_R$)")
ax.set_xlabel("Log$_{10}$ Galaxy Distance (pc)")
#ax.set_title("Galaxy Colour vs Distance")

fig.savefig(dir_path + '/Colour-vs-Distance.png', dpi=200)

plt.clf()

#this plots the galaxy radii vs the mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterRadii), s=3)
ax.set_ylabel("log$_{10}$ Galaxy Radius (pc)")
ax.set_xlabel("log$_{10}$ Galaxy Mass (kg)")
# x = arange(min(log10(clusterMasses)), max(log10(clusterMasses)), 0.1)
# z,cov = polyfit(log10(clusterMasses), log10(clusterRadii), 1, cov=True)     #this finds the linear fit for the data
# p = poly1d(z)
# gradUnc, intUnc = sqrt(diag(cov))
# upper = (z[0] + gradUnc) * x + (z[1] - intUnc)
# lower = (z[0] - gradUnc) * x + (z[1] + intUnc)
# plt.plot(log10(clusterMasses),p(log10(clusterMasses)),"r--", linewidth=0.5)
# plt.fill_between(x, lower, upper, color='r', alpha=0.2)

# print(f"log10(R)=({z[0]:0.5f} \pm {round(gradUnc,5)}) M + ({z[1]:0.3f} \pm {round(intUnc, 2)}) km/s \nR^2 = {r2_score(log10(clusterRadii),p(log10(clusterMasses))):0.3f}")
    
plt.errorbar(log10(clusterMasses), log10(clusterRadii), xerr=logMassUnc, yerr=radUnc, fmt=',', linewidth=0.5)
#ax.set_title("Galaxy Radius vs Mass")

fig.savefig(dir_path + '/Radius vs Mass.png', dpi=300)
    
plt.clf()

#this plots the cluster rotational velocity vs mass
fig, ax = plt.subplots()
plt.scatter(log10(clusterMasses), log10(clusterMaxVel), s=4)
ax.set_ylabel("log$_{10}$ Galaxy Rot Velocity (km s$^{-1}$)")
ax.set_xlabel("log$_{10}$ Galaxy Mass (kg)")

x = arange(min(log10(clusterMasses)), max(log10(clusterMasses)), 0.1)
z,cov = polyfit(log10(clusterMasses), log10(clusterMaxVel), 1, cov=True)     #this finds the linear fit for the data
p = poly1d(z)
gradUnc, intUnc = sqrt(diag(cov))
upper = (z[0] + gradUnc) * x + (z[1] - intUnc)
lower = (z[0] - gradUnc) * x + (z[1] + intUnc)
plt.plot(x, p(x),"r--", linewidth=0.5)
plt.fill_between(x, lower, upper, color='r', alpha=0.2)

print(f"v_r=({z[0]:0.5f} \pm {round(gradUnc,5)}) M + ({z[1]:0.3f} \pm {round(intUnc, 2)}) km/s \nR^2 = {r2_score(log10(clusterMaxVel),p(log10(clusterMasses))):0.3f}")
MassVelError = []
for vel in clusterMaxVel:
    MassVelError.append(0.1 / vel)
plt.errorbar(log10(clusterMasses), log10(clusterMaxVel), xerr=logMassUnc, yerr=MassVelError, fmt=',', linewidth=0.5)
#ax.set_title("Galaxy Rot Velocity vs Mass")
fig.savefig(dir_path + '/Rot Velocity vs Mass.png', dpi=300)

plt.clf()

#this plots the size of a galaxy vs its distance
fig, ax = plt.subplots()
plt.scatter(log10((galaxData['RadialVelocity'] - 1.643) / -0.00263), log10(((galaxData['RadialVelocity'] - 1.643) / -0.00263) * tan(galaxData['Size'])), s=0.1)
ax.set_ylabel("log$_{10}$ Galaxy Size (pc)")
ax.set_xlabel("log$_{10}$ Galaxy Distance (pc)")
#ax.set_title("Galaxy Size vs Distance")

fig.savefig(dir_path + '/Size vs Distance.png', dpi=200)

plt.clf()