# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:42:33 2022

@author: Ryan White     s4499039

This program searches through the *total* fuzzy data to identify clusters of fuzzy objects within the sky. 

To do this, it checks whether fuzzy objects are similar in position in the sky (equatorial and polar angle similarity), as well as if their radial velocities
are similar. The 'closeness' that fuzzy objects have to have is related to their radial velocity, by the function

closeness = 10 * exp(-1/500 * abs(velocity)) + 2

which gives a needed closeness of <=12 degrees for extremely close clusters, and around 3 degrees for clusters with 1000km/s radial velocity. A new function
called closeness() was defined to account for the equatorial angle warping close to polar angles of 0 and 180 (north and south pole respectively).


Cluster/Group Naming Convention:
     - The first character is either 'A' for approaching, or 'R' for receding
     - The second character refers to whether it is a 'C' for cluster, or 'G' for group. For this assignment, a cluster has N >= 20 members. 
         A Group will have 3 <= N < 20 members. This distinction is arbitrary and can be changed via the 'groupthreshold' and 'clusterthreshold' variables. 
     - The first string of characters refers to the objects equatorial angle. (rounded to nearest degree of the central member)
     - The second string of characters refers to the objects polar angle. (rounded to nearest degree of the central member)
     - e.g. RG094-003 would be a receding group, at equatorial angle 94 degrees and polar angle 3 degrees

The output is a file "cluster data.txt" which has the name, number of members, approximate radial velocity of each group/cluster, and the population names of the cluster in the format:
    
    Name    No. of Members  Radial Velocity     Names of Cluster Members
"""

from numpy import *
from math import isclose
import os 
import statistics
import pandas as pd


def closeness(polar, velocity, warp):
    ''' If warp = True, this function outputs an equatorial angle domain such that two galaxies can be considered clusters, given the polar angle and their velocity. 
    If warp = False, this function just gives a polar angle range such that two galaxies can be considered part of the same cluster. 
    '''
    if warp:        #warp is true when looking at equatorial angle
        if 45 <= polar <= 135:
            domain = (10*exp(-1/500 * abs(velocity)) + 2)   #assumed minimal warping for +/-45 degrees either side of the equator
        else:
            domain = (10*exp(-1/500 * abs(veloc)) + 2) / (sin(abs(polar) * pi / 180))   #this accounts for the equatorial angle warping near the poles in the sky
    else:
        domain = (10*exp(-1/500 * abs(veloc)) + 2)      #this is used for the polar angle closeness 
    return domain
    

dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the Fuzzy Cluster Identifier.py directory
totalfuzzy = open(dir_path + "/total fuzzy data.txt", "r")      
fuzzydata = totalfuzzy.readlines()[1:]

for filename in os.listdir(dir_path + "\\Fuzzy Clusters\\"):
    os.remove(dir_path + f'\Fuzzy Clusters\\' + filename)           #removes all files in Fuzzy Clusters folder to prevent double ups on subsequent runs of code

clusters = []       #initializes essential lists for the program
checklist = []

groupthreshold = 3          #define the sizes for groups/clusters
clusterthreshold = 20
i = 0       #main while loop iterating variable
missing = 0

while i <= len(fuzzydata)-1:
    [name, equat, polar, bluef, greenf, redf, size, veloc, dist, location] = fuzzydata[i].split()         #imports a row of data
    
    #following statements clean up the data a bit     
    bluef, greenf, redf, veloc= float(bluef), float(greenf), float(redf), float(veloc)
    equat, polar, dist = float(equat), float(polar), float(dist)
    
    members = []
    workingcluster = []         #initializes the 'current' cluster to focus on, to which individual fuzzy objects are added
    workingdata = []
    j = 0       #secondary while loop iterating variable
    
    while j >= 0:
        if j == 0:      #if this is the first fuzzy to look at
            addendum =  str([name, equat, polar, size, veloc, dist, location]) + "\n"
            workingcluster.append(addendum)
            members.append(name)
            workingdata.append([name, equat, polar, bluef, greenf, redf, size, veloc, dist, location])
            checklist.append(1)
            j += 1
        elif i + j >= len(fuzzydata):       #if the current working cluster index is >= the total data length
            #this code adds the working cluster to the list of clusters, since it has reached it's maximum length
            if len(workingcluster) >= groupthreshold:       #checks if working cluster qualifies as a group
                if sign(veloc) < 0:         #if moving away, cluster begins with "Receding"
                    clustername = "R"
                else:                       #if moving towards us, cluster begins with "Approaching
                    clustername = "A"
                if len(workingcluster) >= clusterthreshold:         #checks if working cluster qualifies as a cluster
                    clustername += "C"          #names it C for cluster
                else:
                    clustername += "G"          #names it G for group
                    
                #following code searches for the middle component of the cluster for naming reasons
                middle = statistics.median(range(1, len(workingcluster)))
                [nameM, equatM, polarM, bluefM, greenfM, redfM, sizeM, velocM, distM, locationM] = fuzzydata[i+round(middle)].split()
                equatM, polarM, distM = float(equatM), float(polarM), float(distM)
                clustername += str('-X' + '%03d' % round(equatM)) + "-Y" + str('%03d' % round(polarM)) + '-N' + str('%03d' % len(workingcluster))      #adds 3 digit, rounded equat and polar angles to clustername
                
                i += j   
                thiscluster = pd.DataFrame(workingdata, columns=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Size', 'RadialVelocity', 'Distance', 'Location'])
                if os.path.exists(dir_path+f'/Fuzzy Clusters/{clustername}.txt'):   #checks if there's another cluster with the same name, adds a '(1)' to clustername if so to differentiate duplicates
                    clusters.append([clustername + '(1)', len(workingcluster), veloc])  #adds the working cluster to the list of clusters
                    thiscluster.to_csv(dir_path+f'/Fuzzy Clusters/{clustername}(1).txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername
                else:
                    clusters.append([clustername, len(workingcluster), veloc])  #adds the working cluster to the list of clusters
                    thiscluster.to_csv(dir_path+f'/Fuzzy Clusters/{clustername}.txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername
                j = -1
            #if the working cluster doesn't qualify as a group or cluster, the else statement then just ignores the working cluster
            else:
                i += j
                missing += len(workingcluster)          #adds the number of fuzzies that aren't part of a group to a counter
                j = -1
            
        else:
            #following code reads the data for a fuzzy object i+j =/= i, and cleans up the data a bit
            [name2, equat2, polar2, bluef2, greenf2, redf2, size2, veloc2, dist2, location2] = fuzzydata[i+j].split()
            veloc2 = float(veloc2)
            equat2, polar2, dist2 = float(equat2), float(polar2), float(dist2)
            
            #following if statement checks if the (i+j)th fuzzy object is close to the ith fuzzy in the sky, AND has a very similar radial velocity
            if (isclose(equat, equat2, abs_tol=closeness(polar, veloc, True)) and isclose(polar, polar2, abs_tol=closeness(polar, veloc, False))) and isclose(veloc, veloc2, rel_tol=0.01):
                addendum = str([name2, equat2, polar2, size2, veloc2, dist2, location2]) + "\n"
                members.append(name2)
                workingcluster.append(addendum)     #if the two objects are close, it adds the (i+j)th object to the working cluster
                checklist.append(1)         #checks the object off so that it isn't read again
                workingdata.append([name2, equat2, polar2, bluef2, greenf2, redf2, size2, veloc2, dist2, location2])
                j += 1
            else:
                #this code reads as the same above
                if len(workingcluster) >= groupthreshold:
                    if sign(veloc) < 0:
                        clustername = "R"
                    else:
                        clustername = "A"
                    if len(workingcluster) >= clusterthreshold:
                        clustername += "C"
                    else:
                        clustername += "G"
                    middle = statistics.median(range(1, len(workingcluster)))
                    [nameM, equatM, polarM, bluefM, greenfM, redfM, sizeM, velocM, distM, locationM] = fuzzydata[i+round(middle)].split()
                    equatM, polarM, distM = float(equatM), float(polarM), float(distM)
                    clustername += str('-X' + '%03d' % round(equatM)) + "-Y" + str('%03d' % round(polarM)) + '-N' + str('%03d' % len(workingcluster))
                    i += j
                    
                    thiscluster = pd.DataFrame(workingdata, columns=['Name', 'Equatorial', 'Polar', 'BlueFlux', 'GreenFlux', 'RedFlux', 'Size', 'RadialVelocity', 'Distance', 'Location'])
                    if os.path.exists(dir_path+f'/Fuzzy Clusters/{clustername}.txt'):  #checks if there's another cluster with the same name, adds a '(1)' to clustername if so to differentiate duplicates
                        clusters.append([clustername + '(1)', len(workingcluster), veloc])  #adds the working cluster to the list of clusters
                        thiscluster.to_csv(dir_path+f'/Fuzzy Clusters/{clustername}(1).txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername
                    else:
                        clusters.append([clustername, len(workingcluster), veloc])  #adds the working cluster to the list of clusters
                        thiscluster.to_csv(dir_path+f'/Fuzzy Clusters/{clustername}.txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername
                    j = -1
                else:
                    i += j      
                    missing += len(workingcluster)
                    j = -1
            
totalfuzzy.close()


#since the above code runs into issues due to the way the data is organised, the following code checks for duplicate clusters and merges them

k = 1
while k <= 3:       #while loop here so that multiple passes of the following error-fixing code are run. k is just the number of passes to use
    if k == 2:
        print("Second pass:")
    if k == 3:
        print("Third pass:")
    for idx, row in enumerate(clusters):            #"for cluster 1"
        [name, n, veloc] = row          #reads data of cluster and then cleans it up a bit
        prefix, equat, polar, population = name.split("-")
        if population[-3:] == '(1)':        #checks if this cluster is a 'duplicate' cluster, and removes the (1) if so
            population = population[:-3]
        equat = equat[-3:]; polar = polar[-3:]; population = population[-3:]
        numequat, numpolar, population = float(equat), float(polar), float(population)
        for index, line in enumerate(clusters):         #"for cluster 2"
            [name2, n2, veloc2] = line
            prefix2, equat2, polar2, population2 = name2.split("-")
            if population2[-3:] == '(1)':
                population2 = population2[:-3]
            equat2 = equat2[-3:]; polar2 = polar2[-3:]; population2 = population2[-3:]
            numequat2, numpolar2, population2 = float(equat2), float(polar2), float(population2)
            #checks if cluster 2 is eerily close to cluster 1 in terms of characteristics. If so, they get merged and their name changed if they now qualify as a cluster
            if (idx != index) and ((numequat2 - 2 <= numequat <= numequat2 + 2) and (numpolar2 - 2 <= numpolar <= numpolar2 + 2) and (veloc2 - 1 <= veloc <= veloc2 + 1)):        #we only care about different clusters
                
                print("MERGING:" + str(clusters[idx]))
                
                if sign(veloc) < 0:
                    clustername = "R"
                else:
                    clustername = "A"
                if n + n2 >= clusterthreshold:
                    clustername += "C"
                else:
                    clustername += "G"
                new_name = clustername + "-X" + equat + "-Y" + polar + "-N" + str('%03d' % int(population + population2))
                clusters[idx] = [new_name, n+n2, veloc] 
                print("DELETING:" + str(clusters[index]))
                print("NEW CLUSTER: ", new_name)
                population += population2   #this and the next line change the main cluster's population to reflect the merger
                n += n2
                with open(dir_path + f'\Fuzzy Clusters\{name}.txt', 'a') as outfile:    #opens main cluster's .txt file
                        with open(dir_path + f'\Fuzzy Clusters\{name2}.txt') as infile: #opens secondary cluster's .txt file
                            for count, line in enumerate(infile):   #for each line in the secondary cluster,
                                if count == 0:
                                    pass
                                else:
                                    outfile.write(line) #adds the secondary cluster's line to the main cluster file
                os.remove(dir_path + f'\Fuzzy Clusters\{name2}.txt')        #removes the 'deleted' cluster which has been merged into the main cluster
                os.rename(dir_path + f'\Fuzzy Clusters\{name}.txt', dir_path + f'\Fuzzy Clusters\{new_name}.txt') #renames the main cluster to reflect the merger
                del clusters[index]     #deletes cluster 2 from the cluster list so that it isn't checked again
                
            else:
                pass
    k += 1

#following code creates a new file "fuzzy clusters.txt" (if it isn't there already) and adds each cluster to it separated by a new line 

allclusters = pd.DataFrame(clusters, columns=['Name', 'No.Members', 'RadialVelocity'])
allclusters.to_csv(dir_path+f'/fuzzy clusters.txt', index=None, sep=' ')    #writes currentcluster to file defined by clustername