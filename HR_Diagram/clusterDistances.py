"""
Script to apply previously written HR.py to star cluster data.
We will read star clusters from directory, then write to csv.
    If star not already in csv, then set distance to 1.
Read off csv distance and plot against baseline.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import logging

import HR as HR

__author__ = "Ciaran Komarakul-Greene"

logging.basicConfig(level=logging.INFO)
dir_path = os.path.dirname(os.path.realpath(__file__))

# get list of all clusters
cluster_dir = dir_path + "/../Data/Camera Images/Star Clusters/"

starClusters = pd.read_csv(dir_path + '/clusterDistances.csv')
# csv of format: (example)
#    ClusterName, Distances
#    X18.3-Y15.4-N101, 132

# make sure all clusters accounted for
# build up list of clusters that are not in the csv already
newClusters = []
for clusterFilename in os.listdir(cluster_dir):
    clusterName = clusterFilename.strip('.txt')
    if not clusterName in list(starClusters['ClusterName']):
        newClusters.append(clusterName)

# if there are new clusters, append to csv file with a rudimentary guess
# at distance; I will match top to top from baseline below
dir_path = os.path.dirname(os.path.realpath(__file__))
# load all stars
calStars = pd.read_csv(dir_path + '/../Data/Camera Images/calibrated star data.txt',
    delimiter=' ')
# filter stars according to parallax
filtStars = calStars.loc[calStars['Parallax'] > HR.PARALLAX_LIMIT] 
filtStarGreen = np.array(filtStars['GreenFlux']) # green flux
filtStarDists = np.array(filtStars['Distance']) # distances
baseline_minMG = min(HR.appToAbs(HR.fluxToApp(filtStarGreen), filtStarDists))

# here we append dataframe with list of all new clusters
if newClusters != []:
    newClusters_df = pd.DataFrame(zip(newClusters, np.ones(len(newClusters))), 
        # create dataframe to concatenate
        columns=['ClusterName', 'Distances', 'Distance Lower Bound', 'Distance Upper Bound']) 
    starClusters = pd.concat((starClusters, newClusters_df), ignore_index=True)

# for each star cluster, plot HR diagram, then save in folder
for clusterName in list(starClusters['ClusterName']):
    # grab data
    clusterData = pd.read_csv(cluster_dir + '/' + clusterName + '.txt', delimiter=' ')
    clusterRed = np.array(clusterData['RedFlux'])
    clusterBlu = np.array(clusterData['BlueFlux'])
    clusterGre = np.array(clusterData['GreenFlux'])
    clusterPer = np.array(clusterData['Periodicity'])

    # here we implement rudimentary estimate 
    if clusterName in newClusters:
        # get cluster max M_green
        cluster_minMG = min(HR.fluxToApp(clusterGre))
        # compare with baseline min and change to distance estimate
        clusterDist = HR.distModToDist(cluster_minMG - baseline_minMG)
        # save to dataframe
        # from HR.ipynb, we have relative uncertainty estimates; constants given in HR.py
        starClusters.loc[starClusters['ClusterName'] == clusterName] = \
            [clusterName, clusterDist, HR.UNCERTAINTY_LOWER*clusterDist, HR.UNCERTAINTY_UPPER*clusterDist]
    else:
        # distance already exists, grab distance
        clusterDist = starClusters.loc[starClusters['ClusterName'] == clusterName].iloc[0,1]

    fig, ax = HR.plotBaseline()
    # compute apparent magnitudes, unpack as arguments, then plot on same figure
    HR.plotHR(*list(map(lambda x: HR.appToAbs(HR.fluxToApp(x), clusterDist),
                [clusterRed, clusterGre, clusterBlu])), clusterPer, ax,col='red',
                lab=clusterName, scale=10)
    ax.set_title(f'Cluster {clusterName} HR Diagram, distance = {clusterDist:.1f} pc')
    ax.legend()
    ax.grid()
    dpi_scale = 2
    fig.savefig(dir_path + f'/ClusterFigures/HR_{clusterName}.png', dpi=dpi_scale*100)
    plt.close(fig)

# save new copy of starClusters for future reference
starClusters.to_csv(dir_path + '/clusterDistances.csv',index=False)

logging.info('done')
