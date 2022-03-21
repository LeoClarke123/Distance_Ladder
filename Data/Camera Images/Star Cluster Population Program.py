"""
Created on Thu Mar  3 15:42:33 2022

@author: Ryan White     s4499039

This program sorts through the "total point-like data.txt" file in relation with the "star clusters.txt" file to identify the stars which are within each star cluster.
This program outputs a .txt document for each star cluster with the data for all of the stars within the cluster (in the same format as "total point-like data.txt").

Of course, this program requires those two files in the same working directory. This program also names the star clusters. Currently, the naming convention is:
    X{coord}-Y{coord}-N{population}.txt
"""

from numpy import *
from astropy.table import Table
import os 
import pandas as pd

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

dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the working directory

clusters = pd.read_csv(dir_path+'/star clusters.txt', delimiter=' ')
stars = pd.read_csv(dir_path+'/total point-like data.txt', delimiter=' ')

cluster_xlocs = clusters['x']       #x and y coords of the center of each cluster
cluster_ylocs = clusters['y']
cluster_x_extents = clusters['horizontal']  #horizontal and vertical extent of the cluster (purely along each axis)
cluster_y_extents = clusters['vertical']
cluster_angles = clusters['theta']      #angle of rotation of cluster relative to positive x-direction

star_equat = stars['Equatorial']
star_polar = stars['Polar']


for cluster in range(len(clusters)):
    data = []
    h = 0
    p = 0
    xcenter = cluster_xlocs[cluster]; ycenter = cluster_ylocs[cluster]; theta = -cluster_angles[cluster] #angle in deg, and negative due to the flipped y-axis (polar coords)
    major = 2 * cluster_x_extents[cluster] / cos(theta*pi/180)  #calculates length of major axis from the x-extent and the angle (converted to radians)
    minor = 2 * cluster_y_extents[cluster] / sin(pi/2 - theta*pi/180)   #calculates length of minor axis from the y-extent and the angle (converted to radians)
    for star in range(len(stars)):
        if (xcenter == 202.5 and ycenter == 140):     #anamolous case for overlapping clusters
            if in_ellipse(211.4, 143.3, 8.61, 2.27, -28, star_equat[star], star_polar[star]): #checks if star is in visually coincident cluster
                pass            #if so, does nothing
            elif in_ellipse(xcenter, ycenter, major, minor, theta, star_equat[star], star_polar[star]):     #checks if star is within parent cluster (the one we care about)
                data.append(stars.loc[star,:])  #adds the star to parent cluster
        elif in_ellipse(xcenter, ycenter, major, minor, theta, star_equat[star], star_polar[star]):
            #add entire row of star data from "stars" here to the current cluster
            data.append(stars.loc[star,:])
        else:
            pass
    #next line converts data list to Pandas data frame, with those column names which match the 'total point-like data.txt' file
    currentcluster = pd.DataFrame(data, columns=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Parallax', 'RadialVelocity', 'Distance', 'Periodicity', 'Location'])
    clustername = 'X' + str(xcenter) + '-Y' + str(ycenter) + '-N' + str(len(data))  #according to naming convention in preamble
    currentcluster.to_csv(dir_path+f'/Star Clusters/{clustername}.txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername
    




