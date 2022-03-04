"""
Created on Thu Mar  3 15:42:33 2022

@author: Ryan White     s4499039

This program searches through the *total* point-like data to identify clusters of point-like objects within the sky. 

To do this, it checks whether point-like objects are similar in position in the sky (equatorial and polar angle similarity), as well as if their radial velocities
are similar. 

Cluster/Group Naming Convention:
     - The first character is either 'A' for approaching, or 'R' for receding
     - The second character refers to whether it is a 'C' for cluster, or 'G' for group. For this assignment, a cluster has N >= 20 members. 
         A Group will have 3 <= N < 20 members. This distinction is arbitrary and can be changed via the 'groupthreshold' and 'clusterthreshold' variables. 
     - The first string of characters refers to the objects equatorial angle. (rounded to nearest degree of the central member)
     - The second string of characters refers to the objects polar angle. (rounded to nearest degree of the central member)
     - e.g. RG094-003 would be a receding group, at equatorial angle 94 degrees and polar angle 3 degrees

The output is a file "star cluster data.txt" which has the name, number of members, and approximate radial velocity of each group/cluster.
"""

from numpy import *
from math import isclose
import os 
import statistics


dir_path = os.path.dirname(os.path.realpath(__file__))      #finds the Cluster Identifier.py directory
totalpoints = open(dir_path + "/total point-like data.txt", "r")      
pointdata = totalpoints.readlines()[1:]

clusters = []       #initializes essential lists for the program
checklist = []

groupthreshold = 3          #define the sizes for groups/clusters
clusterthreshold = 20
i = 0       #main while loop iterating variable
missing = 0


while i <= len(pointdata)-1:
    [name, equat, polar, bluef, greenf, redf, parallax, veloc, distance, periodicity, location] = pointdata[i].split()         #imports a row of data
    
    #following statements clean up the data a bit
    name, location = name.replace('[', ''), location.replace(']', '')       
    bluef, greenf, redf, veloc= float(bluef), float(greenf), float(redf), float(veloc)
    equat, polar = float(equat), float(polar)
    
    workingcluster = []         #initializes the 'current' cluster to focus on, to which individual fuzzy objects are added
    j = 0       #secondary while loop iterating variable
    
    while j >= 0:
        if j == 0:      #if this is the first fuzzy to look at
            addendum =  str([name, equat, polar, parallax, veloc, distance, periodicity, location]) + "\n"
            workingcluster.append(addendum)
            checklist.append(1)
            j += 1
        elif i + j >= len(pointdata):       #if the current working cluster index is >= the total data length
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
                [nameM, equatM, polarM, bluefM, greenfM, redfM, sizeM, velocM, locationM] = pointdata[i+round(middle)].split()
                equatM, polarM = float(equatM), float(polarM)
                clustername += str('%03d' % round(equatM)) + "-" + str('%03d' % round(polarM))      #adds 3 digit, rounded equat and polar angles to clustername
                
                i += j
                clusters.append([clustername, len(workingcluster), veloc])      #adds the working cluster to the list of clusters
                j = -1
            #if the working cluster doesn't qualify as a group or cluster, the else statement then just ignores the working cluster
            else:
                i += j
                missing += len(workingcluster)          #adds the number of fuzzies that aren't part of a group to a counter
                j = -1
            
        else:
            #following code reads the data for a fuzzy object i+j =/= i, and cleans up the data a bit
            [name2, equat2, polar2, bluef2, greenf2, redf2, parallax2, veloc2, distance2, periodicity2, location2] = pointdata[i+j].split()
            veloc2 = float(veloc2)
            equat2, polar2 = float(equat2), float(polar2)
            
            #following if statement checks if the (i+j)th fuzzy object is close to the ith fuzzy in the sky, AND has a very similar radial velocity
            if (isclose(equat, equat2, abs_tol=15) and isclose(polar, polar2, abs_tol=15)) and isclose(veloc, veloc2, rel_tol=0.2):
                addendum = str([name2, equat2, polar2, parallax2, veloc2, distance2, periodicity2, location2]) + "\n"
                workingcluster.append(addendum)     #if the two objects are close, it adds the (i+j)th object to the working cluster
                checklist.append(1)         #checks the object off so that it isn't read again
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
                    [nameM, equatM, polarM, bluefM, greenfM, redfM, parallaxM, velocM, distance2, periodicity2, locationM] = pointdata[i+round(middle)].split()
                    equatM, polarM = float(equatM), float(polarM)
                    clustername += str('%03d' % round(equatM)) + "-" + str('%03d' % round(polarM))
                    i += j
                    clusters.append([clustername, len(workingcluster), veloc])
                    j = -1
                else:
                    i += j      
                    missing += len(workingcluster)
                    j = -1
            
totalpoints.close()

#since the above code runs into issues due to the way the data is organised, the following code checks for duplicate clusters and merges them

k = 1
while k <= 3:       #while loop here so that multiple passes of the following error-fixing code are run. k is just the number of passes to use
    if k == 2:
        print("Second pass:")
    elif k == 3:
        print("Third pass:")
    for idx, row in enumerate(clusters):            #"for cluster 1"
        [name, n, veloc] = row          #reads data of cluster and then cleans it up a bit
        prefix, polar = name.split("-")
        equat = prefix[-3:]
        numequat, numpolar = float(equat), float(polar)
        for index, line in enumerate(clusters):         #"for cluster 2"
            if idx != index:        #we only care about different clusters
                [name2, n2, veloc2] = line
                prefix2, polar2 = name2.split("-")
                equat2 = prefix2[-3:]
                numequat2, numpolar2 = float(equat2), float(polar2)
                #checks if cluster 2 is eerily close to cluster 1 in terms of characteristics. If so, they get merged and their name changed if they now qualify as a cluster
                if (numequat2 - 5 <= numequat <= numequat2 + 5) and (numpolar2 - 5 <= numpolar <= numpolar2 + 5) and (veloc2 - 0.2 <= veloc <= veloc2 + 0.2):
                    print("MERGING:" + str(clusters[idx]))
                    if sign(veloc) < 0:
                        clustername = "R"
                    else:
                        clustername = "A"
                    if n + n2 >= clusterthreshold:
                        clustername += "C"
                    else:
                        clustername += "G"
                    name = clustername + equat + "-" + polar
                    clusters[idx] = [name, n+n2, veloc] 
                    print("DELETING:" + str(clusters[index]))
                    del clusters[index]     #deletes cluster 2 from the cluster list so that it isn't checked again
    k += 1

#following codes creates a new file "cluster data.txt" (if it isn't there already) and adds each cluster to it separated by a new line
totalclusters = open(dir_path+"/WIP star cluster data.txt", "w")
totalclusters.write("Name   No. of Members  Radial Velocity \n")
for cluster in clusters:
    totalclusters.write(str(cluster) + "\n")
