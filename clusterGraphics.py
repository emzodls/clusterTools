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

from fractions import Fraction
import colorsys
from itertools import chain,count
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
'''
color generation from:
http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors
'''
def zenos_dichotomy():
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in count():
        yield Fraction(1,2**k)

def getfracs():
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield 0
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5) # can be used for the v in hsv to map linear values 0..1 to something that looks equidistant

def genhsv(h):
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range

genrgb = lambda x: colorsys.hsv_to_rgb(*x)

flatten = chain.from_iterable

def _get_colors_Janus(num_colors):
    fracGen = getfracs()
    fracs = [fracGen.next() for i in xrange((num_colors+1)/2)]
    rgbs = map(genrgb,flatten(map(genhsv,fracs)))
    return rgbs[:num_colors]

def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


def generateClusterCompGraphic(cluster1,cluster2,pairs,outname):
    noPair = tuple(map(lambda x: x/255.,(211,211,211)))
    distDict = {}
    colorDict = {}
    colorPalette = [tuple(float(value) for value in values) for values in _get_colors_Janus(len(pairs))]

    cluster1Name = "%s:%i-%i" % (cluster1[0].species,cluster1.location[0],cluster1.location[1])
    cluster2Name = "%s:%i-%i" % (cluster2[0].species,cluster2.location[0],cluster2.location[1])

    # Generate color dictionary
    maxLen = 0
    proteinHits = set()
    pairs.reverse()
    for idx,(dist,prot1,prot2) in enumerate(pairs):
        distDict[prot1.hitName] = dist
        distDict[prot2.hitName] = dist
        colorDict[prot1.hitName] = colorPalette[idx]
        colorDict[prot2.hitName] = colorPalette[idx]
        proteinHits.update([prot1,prot2])

    # Draw the cluster
    gd_diagram = GenomeDiagram.Diagram(outname)
    featureHandles = {}
    for idx,cluster in enumerate([cluster1,cluster2]):
        offset = cluster.location[0]
        maxLen = max(maxLen,cluster.location[1]-offset)
        clusterName = "%s:%i-%i" % (cluster[0].species,cluster.location[0],cluster.location[1])
        gd_track_for_features = gd_diagram.new_track(3-2*idx,name = clusterName,
                            start=0, end=cluster.location[1]-offset,
                            scale_ticks=0,scale=0)
        assert clusterName not in featureHandles
        featureHandles[clusterName] = gd_track_for_features.new_set()

    for dist,prot1,prot2 in pairs:
        color = colors.linearlyInterpolatedColor(colors.firebrick,colors.white,  0, 1, dist)
        border = colors.lightgrey

        coord1,direction1 = prot1.location
        coord2,direction2 = prot2.location
        offset1 = cluster1.location[0]
        offset2 = cluster2.location[0]
        coord1 = (x - offset1 for x in coord1)
        coord2 = (x - offset2 for x in coord2)
        F_x = featureHandles[cluster1Name].add_feature(SeqFeature(FeatureLocation(*coord1),strand=0),color=color,border=border)
        F_y = featureHandles[cluster2Name].add_feature(SeqFeature(FeatureLocation(*coord2),strand=0),color=color,border=border)

        gd_diagram.cross_track_links.append(CrossLink(F_x,F_y,color,border))



    for name,cluster in zip([cluster1Name,cluster2Name],[cluster1,cluster2]):
        offset = cluster.location[0]
        for protein in cluster:
            coord,direction  = protein.location
            coord = (x-offset for x in coord)
            if direction == '+':
                strand = +1
            else:
                strand = -1
            feature = SeqFeature(FeatureLocation(*coord),strand=strand)
            featureHandles[name].add_feature(feature,sigil="BIGARROW", color=colorDict.get(protein.hitName,noPair),
                                      name = protein.name,label_position="middle",
                                             label=protein in proteinHits,arrowshaft_height=1,
                                             arrowshaft_length = 0.1,label_strand = 1,
                                   label_size=8, label_angle=45)


    tracks = gd_diagram.get_tracks()
    for track in tracks:
        track.height=1
        track.greytrack_fontcolor  = colors.black
        track.greytrack_labels = 1
        track.greytrack = 1
        track.greytrack_fontsize=16
        track.greytrack_font_rotation = 0
        track.axis_labels = 0
    gd_diagram.draw(format="linear", pagesize='A4', fragments=1,
                start=0, end=maxLen)

    gd_diagram.write(outname + ".svg", "SVG")