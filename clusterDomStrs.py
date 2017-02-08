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
import align,aligntools,hungarian
from scipy.optimize import linear_sum_assignment
from itertools import chain,count
from copy import copy


# I had to switch to the linear_sum_assignment module in scipy for the hungarian algorithm. hungarian wasn't
# doing the matching correctly. This is slower though since it's a python implementation.
# May need to port over a different c implementation if speed becomes an issue

def compareDomStrs(cluster1,cluster2,match,mismatch,gap,scale):
    #first index each cluster
    clus1Strs = dict(enumerate(cluster1))
    clus2Strs = dict(enumerate(cluster2))

    clus1Size = len(cluster1)
    clus2Size = len(cluster2)
    scoreMatrix,word2num,num2word = aligntools.buildScoringDictScaled(chain(*(cluster1+cluster2)),match,mismatch,scale)

    alignScores = np.ndarray((clus1Size,clus2Size))

    # score each pairwise alignment of domain strings to populate a alignment scores matrix
    for i,domStr1 in enumerate(cluster1):
        for j,domStr2 in enumerate(cluster2):
            num1 = [word2num[x] for x in domStr1]
            num2 = [word2num[x] for x in domStr2]
            alignScore,a1,a2 = align.align(num1,num2,gap,gap,scoreMatrix,local=True)
            alignScores[i,j] = alignScore

    #prepare scoring matrix for hungarian algorithm: negate scores and pad matrix
    if clus1Size < clus2Size:
        costMatrix = -np.vstack((copy(alignScores),np.zeros((clus2Size-clus1Size,clus2Size))))
    elif clus2Size < clus1Size:
        costMatrix = -np.hstack((copy(alignScores),np.zeros((clus1Size,clus1Size-clus2Size))))
    else:
        costMatrix = -copy(alignScores)

    # apply hungarian algorithm for matching
    pairings = [(x,y) for x,y in zip(*linear_sum_assignment(costMatrix)) if (x<clus1Size) and (y<clus2Size)]
    clusterScore = sum(alignScores[pairing] for pairing in pairings)
    pairStrings = [(alignScores[(x,y)],clus1Strs[x],clus2Strs[y]) for x,y in pairings]
    pairStrings.sort(reverse=True)

    return clusterScore,pairStrings

def parseDomStrsFile(domStrFile):
    domStrDict = {}
    with open(domStrFile,'rb') as domStrHandle:
        for line in domStrHandle:
            if line[0] == '>':
                clusterName = line[1:-1]
                domStrDict[clusterName] = []
            else:
                domStrDict[clusterName].append(line.split())
    return domStrDict