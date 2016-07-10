# Copyright (C) 2016 Emmanuel LC. de los Santos
# University of Warwick
# Warwick Integrative Synthetic Biology Centre
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

'''
    This file is part of clusterTools.

    clusterTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    clusterTools is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with clusterTools.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
from scipy.optimize import linear_sum_assignment

# I had to switch to the linear_sum_assignment module in scipy for the hungarian algorithm. hungarian wasn't
# doing the matching correctly. This is slower though since it's a python implementation.
# May need to port over a different c implementation if speed becomes an issue

def calculateClusterDist(cluster1,cluster2,hitDictID,linearDist = True):
    '''
    Given two clusters annotated with the same hitDict (typically and all-v-all comparison) will estimate the distance
    between the two clusters by calculating the normalized distance between each protein in the clusters. It will use
    this information to create a maximal matching will then scale the matching pairs distance by the percentage of the
    cluster each matching will give. Returns a value between zero and one and the maximal matching of the clusters
    '''

    clus1Size = len(cluster1)
    clus2Size = len(cluster2)

    scoreMatrix =  np.ndarray((clus1Size,clus2Size))

    clus1ProtSize = float(sum(protein.size() for protein in cluster1))
    clus2ProtSize = float(sum(protein.size() for protein in cluster2))

    # populate the score matrix if there are any proteins that are "close together"
    for i,proteinI in enumerate(cluster1):
        for j,proteinJ in enumerate(cluster2):
            scoreMatrix[i,j] = proteinI.calculate_distance(proteinJ,hitDictID,linearDist=linearDist)

    # get the pairings
    pairings = [(x,y) for x,y in zip(*linear_sum_assignment(scoreMatrix)) if (x<clus1Size) and (y<clus2Size)]
    pairScores = [(scoreMatrix[(x,y)],cluster1[x],cluster2[y]) for x,y in pairings]
    pairs = [(x,y,z) for x,y,z in pairScores if x < 1.]
    pairs.sort()
    # scale by the one with less coverage
    clus1cvg = sum(entry[1].size()/clus1ProtSize for entry in pairs)
    clus2cvg = sum(entry[2].size()/clus2ProtSize for entry in pairs)

    if clus1cvg < clus2cvg:
        lessCvgIdx = 1
        lessCvgSize = clus1ProtSize
    else:
        lessCvgIdx = 2
        lessCvgSize = clus2ProtSize
    # scale distances by larger cluster
    # print clus1cvg,clus2cvg
    # print lessCvgIdx,lessCvgSize
    distance = 0
    percentNonHit = 1
    for entry in pairs:
        distance += entry[0]*(entry[lessCvgIdx].size()/lessCvgSize)
        percentNonHit -= (entry[lessCvgIdx].size()/lessCvgSize)
        # print entry[1].hitName,entry[2].hitName,distance, percentNonHit
    percentNonHit = max(0,percentNonHit)
    # print percentNonHit
    distance += percentNonHit

    return distance,pairs

def calculateDistBitScore(cluster1,cluster2,hitDictID):
    '''
    another way to measure distance that will pair up the proteins using matching
    then using the bitscores of the pairs of proteins added up will calculate the distance
    in a similar way to how protein distance is calculated:

    linear: 1 - (sum bit score of paired matches/max bit score of cluster)

    scaling will work this way:

    If there are reciprocal scores average the distance. Otherwise, if hits are on smaller
    cluster, (small clus/large clus)*dist + (1 - small clus/large clus)

    '''
    clus1Size = len(cluster1)
    clus2Size = len(cluster2)

    scoreMatrix =  np.ndarray((clus1Size,clus2Size))

    clus1ProtSize = float(sum(protein.size() for protein in cluster1))
    clus2ProtSize = float(sum(protein.size() for protein in cluster2))

    # populate the score matrix if there are any proteins that are "close together"
    for i,proteinI in enumerate(cluster1):
        for j,proteinJ in enumerate(cluster2):
            scoreMatrix[i,j] = proteinI.calculate_distance(proteinJ,hitDictID,linearDist=True)

    # get the pairings
    pairings = [(x,y) for x,y in zip(*linear_sum_assignment(scoreMatrix)) if (x<clus1Size) and (y<clus2Size)]
    pairScores = [(scoreMatrix[(x,y)],cluster1[x],cluster2[y]) for x,y in pairings]
    pairs = [(x,y,z) for x,y,z in pairScores if x < 1.]
    pairs.sort()
    # figure out which cluster has the hits and calculate coverage

    clus1Flag = sum(1 for protein in cluster1 if protein.hitName in protein.hit_dict[hitDictID].hits) == clus1Size
    clus2Flag = sum(1 for protein in cluster2 if protein.hitName in protein.hit_dict[hitDictID].hits) == clus2Size

    # check that at least one of the clusters has hits to the other cluster
    try:
        assert clus1Flag or clus2Flag
    except AssertionError:
        "print Error: there doesn't seem to be homology information in either of the clusters"

    if clus1Flag and clus2Flag:
        clus1maxBitScore = sum(protein.hit_dict[hitDictID].maxscore for protein in cluster1)
        clus2maxBitScore = sum(protein.hit_dict[hitDictID].maxscore for protein in cluster2)

        clus1cumBitScore = sum(proteinI.hit_dict[hitDictID].get(proteinJ.hitName,0) for dist,proteinI,proteinJ in pairs)
        clus2cumBitScore = sum(proteinJ.hit_dict[hitDictID].get(proteinI.hitName,0) for dist,proteinI,proteinJ in pairs)

        return 0.5*(1 - clus1cumBitScore/clus1maxBitScore) + 0.5*(1-clus2cumBitScore/clus2maxBitScore),pairs
    elif clus1Flag:
        clus1maxBitScore = sum(protein.hit_dict[hitDictID].maxscore for protein in cluster1)
        clus1cumBitScore = sum(proteinI.hit_dict[hitDictID].get(proteinJ.hitName,0) for dist,proteinI,proteinJ in pairs)

        if clus1ProtSize < clus2ProtSize:
            return 1-(clus1cumBitScore/clus1maxBitScore)*(clus1ProtSize/clus2ProtSize),pairs
        else:
            return (1-clus1cumBitScore/clus1maxBitScore),pairs
    else:
        clus2maxBitScore = sum(protein.hit_dict[hitDictID].maxscore for protein in cluster1)
        clus2cumBitScore = sum(proteinI.hit_dict[hitDictID].get(proteinJ.hitName,0) for dist,proteinI,proteinJ in pairs)

        if clus1ProtSize > clus2ProtSize:
            return 1-(clus2cumBitScore/clus2maxBitScore)*(clus2ProtSize/clus1ProtSize),pairs
        else:
            return (1-clus2cumBitScore/clus2maxBitScore),pairs


def getNewick(node, newick, parentdist, leaf_names):
    # from http://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format jfn
    if node.is_leaf():
        return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
