# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 11:29:33 2022

@author: Ryan White
"""

from numpy import *
import os 
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

dir_path = os.path.dirname(os.path.realpath(__file__))


clusterdistances = pd.read_csv(dir_path.replace("Data", "HR_Diagram/clusterDistances.csv"), delimiter=',')

longGradUnc = 0
longIntUnc = 0
shortGradUnc = 0
shortIntUnc = 0
meandifferences = []

for cluster in os.listdir(dir_path + "\\Camera Images\\Star Clusters"):
    stardata = pd.read_csv(dir_path + f"\\Camera Images\\Star Clusters\\{cluster}", delimiter=' ')
    clustername = cluster.replace(".txt", "")
    clusterDist = clusterdistances.loc[clusterdistances["ClusterName"]==clustername].iloc[0,1]
    clusterDistL = clusterdistances.loc[clusterdistances["ClusterName"]==clustername].iloc[0,2]
    clusterDistU = clusterdistances.loc[clusterdistances["ClusterName"]==clustername].iloc[0,3]
    clusterDistUnc = abs(clusterDist - (clusterDistL + clusterDistU)/2)
    starDist = stardata['Distance']
    starPeriod = stardata['Periodicity']
    starGreen = stardata['GreenFlux']
    shortPDists = []; longPDists = []; PDists = []
    differences = []; Ldifferences = []; Sdifferences = []
    for index, distance in enumerate(starDist):
        period = starPeriod[index]
        HRdistance = starDist
        flux = starGreen[index]
        if period != 0 and period > 30:
            periodLDist = sqrt(exp(0.34 * period - 41.3) / flux)
            periodLDistUnc = 2 * periodLDist * sqrt((longGradUnc * period)**2 + longIntUnc**2)
            longPDists.append(periodLDist)
            PDists.append(periodLDist)
            differences.append(periodLDist - clusterDist)
            meandifferences.append(periodLDist - clusterDist)
            Ldifferences.append(periodLDist - clusterDist)
        elif period != 0 and period > 10:
            periodSDist = sqrt(exp(-0.64 * period - 3.3) / flux)
            periodSDistUnc = 2 * periodSDist * sqrt((shortGradUnc * period)**2 + shortIntUnc**2)
            shortPDists.append(periodSDist)
            PDists.append(periodSDist)
            differences.append(periodSDist - clusterDist)
            meandifferences.append(periodSDist - clusterDist)
            Sdifferences.append(periodSDist - clusterDist)
    
    if len(PDists) > 5:
        # meandifferences.append(differences)
        print(clustername, median(Ldifferences), "and", std(Ldifferences))
    # print(mean(differences), std(differences))
    
    
print("Mean diff = ", median(meandifferences), "with SD = ", std(meandifferences))
