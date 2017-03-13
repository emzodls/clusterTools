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

from math import log
from copy import copy,deepcopy
from collections import defaultdict
from sortedcontainers import SortedDict, SortedListWithKey
from operator import itemgetter
from bx.intervals.intersection import IntervalTree, Interval
from Bio import SeqIO
from itertools import count,repeat,groupby
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess,os

class hitDict():
    def __init__(self):
        self.score_dict = dict()
        self.maxscore = 1.
        self.hits = set()

    def add_hit(self,hit_id,score):
        if score > self.maxscore:
            self.maxscore = score
        if hit_id in self.hits:
            if self.score_dict[hit_id] < score:
                self.score_dict[hit_id] = score
        else:
            self.score_dict[hit_id] = score
            self.hits.add(hit_id)
    def get(self,key,default=None):
        return self.score_dict.get(key,default)
    def __repr__(self):
        return repr(self.score_dict)
class Protein():
    # Hit ID will be a tuple (protein ID,internal ID (species_ID_CDS_idx))
    def __init__(self,species_id,protein_id,internalID,idx,location):
        self.species = species_id
        self.name = protein_id
        self.location = location
        self.annotations = dict()
        self.hit_dict = defaultdict(hitDict)
        self.idx = idx
        self.internalID = internalID
        self.hitName = (self.name,self.internalID)
        self.sequence = ''
        self.fastaID = '%s|%i-%i|%s|%s|%s'% (self.species,
                                             self.location[0][0],
                                             self.location[0][1],
                                             self.location[1],
                                             self.internalID,
                                             self.name)
    def __repr__(self):
        return ('Protein(Name=%s,Species=%s,Index=%s,Hit Dicts=%s,Annotations=%s)' % (repr(self.name),
                                                                                      repr(self.species),
                                                                                      repr(self.idx),
                                                                                      repr(self.hit_dict.keys()),
                                                                                      repr(self.annotations.keys())))

    def calculate_distance(self, protein, hit_dict_id,linearDist = True):
        assert (hit_dict_id in self.hit_dict.keys() or hit_dict_id in protein.hit_dict.keys())
        if self.hitName == protein.hitName:
            return 0.
        else:
            if protein.hitName in self.hit_dict[hit_dict_id].hits and \
                            self.hitName not in protein.hit_dict[hit_dict_id].hits:
                if not linearDist:
                    normScore = -log(self.hit_dict[hit_dict_id].get(protein.hitName) / self.hit_dict[hit_dict_id].maxscore)
                    maxDist = -log(1 / self.hit_dict[hit_dict_id].maxscore)
                    if maxDist == 0:
                        return 1
                    else:
                        return normScore / maxDist
                else:
                    distance = 1 - (self.hit_dict[hit_dict_id].get(protein.hitName) / self.hit_dict[hit_dict_id].maxscore)
                    return distance

            elif protein.hitName not in self.hit_dict[hit_dict_id].hits and \
                            self.hitName in protein.hit_dict[hit_dict_id].hits:
                if not linearDist:
                    normScore = -log(protein.hit_dict[hit_dict_id].get(self.hitName) / protein.hit_dict[hit_dict_id].maxscore)
                    maxDist = -log(1 / protein.hit_dict[hit_dict_id].maxscore)
                    if maxDist == 0:
                        return 1
                    else:
                        return normScore / maxDist
                else:
                    distance = 1 - (protein.hit_dict[hit_dict_id].get(self.hitName) / protein.hit_dict[hit_dict_id].maxscore)
                    return distance

            else:
                if len(self.hit_dict[hit_dict_id].hits) == 0  or len(protein.hit_dict[hit_dict_id].hits) == 0:
                    proteinAbundanceFactor = 1
                    targetAbundanceFactor = 1
                else:
                    proteinAbundanceFactor = len(self.hit_dict[hit_dict_id].hits & protein.hit_dict[hit_dict_id].hits) / \
                                     float(len(self.hit_dict[hit_dict_id].hits))
                    targetAbundanceFactor = len(self.hit_dict[hit_dict_id].hits & protein.hit_dict[hit_dict_id].hits) / \
                                    float(len(protein.hit_dict[hit_dict_id].hits))

                # if both the protein abundance factor and the target abundance factor are 0 this means that there
                # aren't any common blast hits between the proteins so the distance is 1
                if proteinAbundanceFactor == 0 and targetAbundanceFactor == 0:
                    return 1.
                elif not linearDist:
                    normScore = self.hit_dict[hit_dict_id].get(protein.hitName, 1) / self.hit_dict[hit_dict_id].maxscore
                    normReciprocalScore = protein.hit_dict[hit_dict_id].get(self.hitName, 1) / protein.hit_dict[hit_dict_id].maxscore
                    maxDist = -log((proteinAbundanceFactor / self.hit_dict[hit_dict_id].maxscore +
                                    targetAbundanceFactor / protein.hit_dict[hit_dict_id].maxscore) / 2)
                    distance = -log(((normScore * proteinAbundanceFactor) + (normReciprocalScore * targetAbundanceFactor)) / 2)
                    if maxDist == 0:
                        return 1
                    else:
                        return distance / maxDist
                else:
                    normScore = protein.hit_dict[hit_dict_id].get(self.hitName, 0) / protein.hit_dict[hit_dict_id].maxscore
                    normReciprocalScore = protein.hit_dict[hit_dict_id].get(self.hitName, 0) / protein.hit_dict[hit_dict_id].maxscore
                    distance = (proteinAbundanceFactor * (1 - normScore)\
                               + targetAbundanceFactor * (1 - normReciprocalScore)) \
                               / (proteinAbundanceFactor+targetAbundanceFactor)
                    return distance

    def size(self):
        return int(calculate_window(self.location[0])/3-1)
    def add_hit(self,hit_dict_id,hit_id,score):
        self.hit_dict[hit_dict_id].add_hit(hit_id,score)
    def extractSequence(self,coordinates):
        try:
            return self.sequence[coordinates[0]-1:coordinates[1]]
        except Exception as exc:
            print('Excepction: %s' % exc.message)
            return ''
    def getAnnotations(self,anotID):
        if anotID in self.annotations.keys():
            return set(x[0] for x in self.annotations[anotID].values())
        else:
            return set()
    def getDomStr(self,anotID,delim):
        if anotID in self.annotations.keys():
            return delim.join(hit[0] for hit in self.annotations[anotID].itervalues())
        else:
            return ''

class Cluster(SortedListWithKey):
    def __init__(self,proteins):
        super(Cluster,self).__init__(proteins,key=lambda x: x.location[0])
        self.location = []
        self.species = None
        if len(self) > 0:
            self.location = [self[0].location[0][0],self[-1].location[0][1]]
            self.species = self[0].species
    def __repr__(self):
        return ('Cluster(Size = %s, Proteins = %s)' % (repr(calculate_window(self.location)),
                                                      repr(', '.join(['(%i,%s)' % (x.idx,x.name) for x in self]))))
    def size(self):
        return calculate_window(self.location)
    def add(self,protein):
        super(Cluster,self).add(protein)
        self.location = [self[0].location[0][0],self[-1].location[0][1]]
    def pop(self,idx=-1):
        super(Cluster,self).pop(idx)
        self.location = [self[0].location[0][0],self[-1].location[0][1]]

def resolve_conflicts(pfam_hit_dict,minDomSize = 9,verbose=False):
    '''
    :param pfam_hit_dict: dictionary of hits for the gene in the following format
    hit start,hit end : int
    hit id : str
    score, model coverage percent : float
    {(hit start,hit end):('hit id',score,model coverage percent)}
    :param minDomSize: int, the minimum window size that will be considered a domain
    :return:
    a sorted dictionary with the position of the hit as the keys and ('hit id',score,model coverage percent)
    '''
    # initialize output
    gene_hits = SortedDict()
    redoFlag = True
    while redoFlag:
        if verbose: print("Sorting through intervals", pfam_hit_dict)
        redoFlag = False
        intervals_scores = [(key,value[1]) for key,value in pfam_hit_dict.items()]
        # sort intervals from pfam hits by score and place the highest score first
        intervals_scores.sort(key=itemgetter(1),reverse=True)
        # initialize intersect tree for quick overlap search
        intersectTree = IntervalTree()
        #add the intervals with the highest scores first
        for (interval,score) in intervals_scores:
            intervalStart = interval[0]
            intervalEnd = interval[1]
            intervalLength = intervalEnd-intervalStart+1
            # if the interval is less than the minimum domain size don't bother
            if intervalLength > minDomSize:
                intersectingIntervals = [(x.start,x.end) for x in intersectTree.find(intervalStart,intervalEnd)]
                overLapFlag = False
                # for every interval that you're adding resolve the overlapping intervals
                while len(intersectingIntervals) > 0 and intervalLength > 1:

                    start,end = intersectingIntervals[0]

                    # interval completely covers existing coverage, break up into two intervals and redo the process
                    if (intervalStart <= start and intervalEnd >= end):
                        if verbose: print("Split Interval", interval,intersectingIntervals, pfam_hit_dict[interval])
                        left_scale = calculate_window((intervalStart,start-1))/intervalLength
                        right_scale = calculate_window((end+1,intervalEnd))/intervalLength
                        pfam_hit_dict[(intervalStart,start-1)] = (pfam_hit_dict[interval][0],
                                                                  pfam_hit_dict[interval][1],
                                                                  pfam_hit_dict[interval][2] * left_scale)
                        pfam_hit_dict[(end+1,intervalEnd)] = (pfam_hit_dict[interval][0],
                                                              pfam_hit_dict[interval][1],
                                                              pfam_hit_dict[interval][2] * right_scale)
                        # delete original hit and iterate
                        del pfam_hit_dict[interval]
                        redoFlag = True
                        break
                    else:
                        #completely in the interval
                        if (intervalStart >= start and intervalEnd <= end):
                            #if completely overlapping then ignore since we already sorted by score
                            overLapFlag = True
                            break
                        #intersection covers the left hand side of the interval
                        elif intervalStart >= start:
                            intervalStart = end + 1
                        #intersection covers the right hand side of the interval
                        elif intervalEnd <= end:
                            intervalEnd = start - 1
                            # recalculate the interval length and see if there are still intersecting intervals
                        intervalLength = intervalEnd-intervalStart+1
                        intersectingIntervals = [(x.start,x.end) for x in intersectTree.find(intervalStart,intervalEnd)]

                if redoFlag:
                    if verbose: print("Exiting For Loop to Reinitialize",pfam_hit_dict)
                    break
                # if loop did not break because of an overlap add the annotation after resolving overlap,
                # check for minimum length after you merge intervals
                elif not overLapFlag and intervalLength > minDomSize:
                    if verbose: print("Adding Hit",(intervalStart,intervalEnd),pfam_hit_dict[interval][0])
                    # scale the hitCoverage based on the reduction this works since interval is a tuple and isn't mutated
                    hitCoverage = pfam_hit_dict[interval][2]*(intervalLength/(interval[1]-interval[0]+1.))
                    gene_hits[(intervalStart,intervalEnd)] = (pfam_hit_dict[interval][0],
                                                              pfam_hit_dict[interval][1],
                                                              hitCoverage)
                    intersectTree.add_interval(Interval(intervalStart,intervalEnd))
    if verbose: print("Merging Hits")
    # Merge Windows Right Next to one another that have the same pFam ID,
    # redoFlag: need to restart the process after a successful merge
    redoFlag = True
    while redoFlag:
        for idx in range(len(gene_hits)-1):
            left_hit = gene_hits.keys()[idx]
            right_hit = gene_hits.keys()[idx+1]
            left_window_size = calculate_window(left_hit)
            right_window_size = calculate_window(right_hit)
            merged_window_size = calculate_window((left_hit[0],right_hit[1]))
            new_coverage = (gene_hits[left_hit][2] + gene_hits[right_hit][2])*\
                           (left_window_size+ right_window_size)/merged_window_size
            # Will merge a hit under the following conditions:
            # 1. Gap between the two hits is less than the minimum domain
            # 2. Cumulative coverage of the two hits is less than 1 (this avoids merging repeats together)
            if right_hit[0]-left_hit[1] < minDomSize and gene_hits[left_hit][0] == gene_hits[right_hit][0] \
                    and new_coverage < 1:
                gene_hits[(left_hit[0],right_hit[1])] = (gene_hits[left_hit][0],
                                                         left_window_size/merged_window_size * gene_hits[left_hit][1] +
                                                         right_window_size/merged_window_size * gene_hits[right_hit][1],
                                                         new_coverage)
                redoFlag = True
                del gene_hits[left_hit]
                del gene_hits[right_hit]
                if verbose: print("Merged", left_hit,right_hit)
                break
        else:
            redoFlag = False
    if verbose: print("Deleting Domains Under Minimum Domain Size")
    # Finally check if any of the domains are less than the minimum domain size
    keysToDelete = [coordinates for coordinates in gene_hits.keys() if calculate_window(coordinates) < minDomSize]
    for key in keysToDelete:
        del gene_hits[key]
        if verbose: print("Deleting",key)
    if verbose: print("Final Annotation", gene_hits)
    return gene_hits

def parse_allVall_blast_file(path,proteinDict,swapQuery=False,evalCutoff=10,scoreCutoff=0,speciesFilter = None):
    if swapQuery:
        queryIdx = 1
        hitIdx = 0
    else:
        queryIdx = 0
        hitIdx = 1
    with open(path) as blast_handle:
        for line in blast_handle:
            line_parse = line.split('\t')
            query_parse = line_parse[queryIdx].split('|')
            print(query_parse)
            species_id = query_parse[0]
            if ((not speciesFilter) or (species_id in speciesFilter)):
                coordinates = [int(x) for x in query_parse[1].split('-')]
                direction = query_parse[2]
                queryIntID = query_parse[3]
                queryProtID = query_parse[4]
                queryIntIdx = int(queryIntID.split('_')[-1])
                protein = proteinDict.setdefault(queryIntID,Protein(species_id,queryProtID,queryIntID,
                                                queryIntIdx,(tuple(coordinates),direction)))
                hit_parse = line_parse[hitIdx].split('|')
                print(hit_parse)
                hit_id = (hit_parse[4],hit_parse[3])
                # assuming blast outfmt 6 output
                score  = float(line_parse[11])
                eval = float(line_parse[10])
                if eval <= evalCutoff and score >= scoreCutoff:
                    protein.add_hit('blast',hit_id,score)
    return proteinDict

def parseBLAST(path,proteinDict,swapQuery=False,evalCutoff=10,scoreCutoff=0,speciesFilter = None):
    if swapQuery:
        queryIdx = 1
        hitIdx = 0
    else:
        queryIdx = 0
        hitIdx = 1
    with open(path) as blast_handle:
        for line in blast_handle:
            try:
                line_parse = line.split('\t')
                query_parse = line_parse[queryIdx].split('|')
                species_id = query_parse[0]
                if ((not speciesFilter) or (species_id in speciesFilter)):
                    coordinates = [int(x) for x in query_parse[1].split('-')]
                    direction = query_parse[2]
                    queryIntID = query_parse[3]
                    queryProtID = query_parse[4]
                    queryIntIdx = int(queryIntID.split('_')[-1])
                    protein = proteinDict.setdefault(queryIntID,Protein(species_id,queryProtID,queryIntID,
                                                    queryIntIdx,(tuple(coordinates),direction)))
                    hit_id = line_parse[hitIdx]
                    score  = float(line_parse[11])
                    eval = float(line_parse[10])
                    if eval <= evalCutoff and score >= scoreCutoff:
                        protein.add_hit('blast',hit_id,score)
            except (ValueError,IndexError):
                pass
    return proteinDict

def calculate_window(coordinates):
    return coordinates[-1]-coordinates[0]+1.
def cluster_proteins(proteins,window_size):
    '''
    Given a list of proteins and a window size will return a dictionary whose key is the species and elements are
    clusters of different proteins in the species
    :param proteins: List of Protein objects, window_size
    :return: dictionary {species:[clusters]}
    '''
    cluster_dict = dict()
    for protein in proteins:
        if protein.species in cluster_dict.keys():
            species_clusters = cluster_dict[protein.species]
            _inCluster  = False
            for species_cluster in species_clusters:
                new_coords = list(copy(species_cluster.location))
                new_coords.extend(protein.location[0])
                new_coords.sort()
                # Changed workflow here to allow proteins to get placed in multiple clusters,
                # if no appropriate cluster is found a new cluster is created for the protein
                if calculate_window(new_coords) <= window_size:
                    species_cluster.add(protein)
                    _inCluster = True
            # if none of the clusters in the dictionary meet the criterion, generate a new cluster
            if not _inCluster:
                species_clusters.append(Cluster([protein]))
        else:
            cluster_dict[protein.species] = [Cluster([protein])]
        # for each cluster in the species check if the protein will fit without going over the frame window
    return cluster_dict

def clusterProteins(proteins,windowSize):
    '''
    Given a list of proteins and a window size will return a dictionary whose key is the species and elements are
    clusters of different proteins in the species. This modifies the old cluster proteins method, sorting the proteins by
    species first based on coordinates, and then sorting them together
    :param proteins: List of Protein objects, window_size
    :return: dictionary {species:[clusters]}
    '''
    proteinsBySpecies = dict()
    cluster_dict = dict()
    # group proteins by species
    for protein in proteins:
        speciesProteins = proteinsBySpecies.setdefault(protein.species,[])
        speciesProteins.append(protein)
    # sort each of the species lists by coordinates
    for proteins in proteinsBySpecies.values():
        proteins.sort(key=lambda x: x.location)
    # Now cluster the proteins
    for species,proteins in proteinsBySpecies.items():
        for protein in proteins:
            clusters = cluster_dict.setdefault(species,[Cluster([protein])])
            for cluster in clusters:
                new_coords = list(copy(cluster.location))
                new_coords.extend(protein.location[0])
                new_coords.sort()
                if calculate_window(new_coords) <= windowSize:
                    cluster.add(protein)
            # add protein as a new cluster in case
            clusters.append(Cluster([protein]))
    return cluster_dict


def parse_hmmscan_domtbl_anot(path,minDomSize,anotID,proteinDict,cutoff_score=25,verbose = False):
    """
    :param path: path to pfam hits that you want parsed

    :return: {protein:gene hit dictionary (output of resolve_conflicts)}
    hit dictionary of pfam hits that can be linked to proteins for their pfam domain annotation
    """
    hit_dict = {}
    with open(path) as hmmfile:
        current_gene_id = ''
        for line in hmmfile:
            if line[0] == '#':
                pass
            else:
                hit = line.split()
                score = float(hit[13])
                if score >= cutoff_score:
                    gene_id = hit[3]
                    #check if you're done parsing the hits of the old gene and annotate the gene
                    if gene_id != current_gene_id and current_gene_id != '':
                        current_gene_parse = current_gene_id.split('|')
                        species_id = current_gene_parse[0]
                        coordinates = [int(x) for x in current_gene_parse[1].split('-')]
                        direction = current_gene_parse[2]
                        geneIntID = current_gene_parse[3]
                        geneProtID = current_gene_parse[4]
                        geneIntIdx = int(geneIntID.split('_')[-1])
                        protein = proteinDict.setdefault(geneIntID,Protein(species_id,geneProtID,geneIntID,geneIntIdx,(tuple(coordinates),direction)))
                        #first resolve overlaps with current gene dict
                        hmmer_hits = resolve_conflicts(hit_dict,minDomSize=minDomSize,verbose = verbose)
                        # check if any hits pass the filter
                        if len(hmmer_hits) > 0:
                            protein.annotations[anotID] = hmmer_hits
                        hit_dict = {}
                    current_gene_id = gene_id
                    hmmer_hit = hit[0]
                    hmm_coverage_start = int(hit[15])
                    hmm_coverage_end = int(hit[16])
                    hmm_length = int(hit[2])
                    hmm_coverage = (hmm_coverage_end-hmm_coverage_start+1.)/hmm_length
                    gene_hmm_start = int(hit[17])
                    gene_hmm_end = int(hit[18])
                    hit_dict[(gene_hmm_start,gene_hmm_end)] = (hmmer_hit,score,hmm_coverage)
    return proteinDict

def parse_hmmsearch_domtbl_anot(path,minDomSize,anotID,proteinDict,eval_cutoff=1e-5,cutoff_score=25,verbose = False,
                                speciesFilter = None,proteinFilter = None):
    """
    :param path: path to pfam hits that you want parsed

    :return: {protein:gene hit dictionary (output of resolve_conflicts)}
    hit dictionary of pfam hits that can be linked to proteins for their pfam domain annotation
    """
    hit_dict = {}
    with open(path) as hmmfile:
        for line in hmmfile:
            if line[0] == '#':
                pass
            else:
                hit = line.split()
                score = float(hit[13])
                gene_id = hit[0]
                evalue = float(hit[12])
                gene_id_parse = gene_id.split('|')
                species_id = gene_id_parse[0]
                if ((not speciesFilter) or (species_id in speciesFilter)) and \
                        ((not proteinFilter) or (gene_id in proteinFilter)):
                    if score >= cutoff_score  and evalue <= eval_cutoff:
                        #check if you're done parsing the hits of the old gene and annotate the gene
                        coordinates = [int(x) for x in gene_id_parse[1].split('-')]
                        direction = gene_id_parse[2]
                        geneIntID = gene_id_parse[3]
                        geneProtID = gene_id_parse[4]
                        geneIntIdx = int(geneIntID.split('_')[-1])
                        protein = proteinDict.setdefault(geneIntID,Protein(species_id,geneProtID,geneIntID,geneIntIdx,(tuple(coordinates),direction)))
                        hit_dict = protein.annotations.get(anotID,{})
                        # check if any hits pass the filter
                        hmmer_hit = hit[3]
                        hmm_coverage_start = int(hit[15])
                        hmm_coverage_end = int(hit[16])
                        hmm_length = int(hit[5])
                        hmm_coverage = (hmm_coverage_end-hmm_coverage_start+1.)/hmm_length
                        gene_hmm_start = int(hit[17])
                        gene_hmm_end = int(hit[18])
                        hit_dict[(gene_hmm_start,gene_hmm_end)] = (hmmer_hit,score,hmm_coverage)
                        #first resolve overlaps with current gene dict
                        hmmer_hits = resolve_conflicts(hit_dict,minDomSize=minDomSize,verbose = verbose)
                        if len(hmmer_hits) > 0:
                            protein.annotations[anotID] = hmmer_hits
    return proteinDict

def parse_phmmer_domtbl_anot(path,minDomSize,anotID,proteinDict,eval_cutoff=1e-5,cutoff_score=25,verbose = False):
    """
    :param path: path to pfam hits that you want parsed

    :return: {protein:gene hit dictionary (output of resolve_conflicts)}
    hit dictionary of pfam hits that can be linked to proteins for their pfam domain annotation
    """
    hit_dict = {}
    with open(path) as hmmfile:
        current_gene_id = ''
        for line in hmmfile:
            if line[0] == '#':
                pass
            else:
                hit = line.split()
                eval = float(hit[12])
                score = float(hit[13])
                if score >= cutoff_score and eval <= eval_cutoff:
                    gene_id = hit[0]
                    #check if you're done parsing the hits of the old gene and annotate the gene
                    if gene_id != current_gene_id and current_gene_id != '':
                        current_gene_parse = current_gene_id.split('|')
                        species_id = current_gene_parse[0]
                        coordinates = [int(x) for x in current_gene_parse[1].split('-')]
                        direction = current_gene_parse[2]
                        geneIntID = current_gene_parse[3]
                        geneProtID = current_gene_parse[4]
                        geneIntIdx = int(geneIntID.split('_')[-1])
                        protein = proteinDict.setdefault(geneIntID,Protein(species_id,geneProtID,geneIntID,geneIntIdx,(tuple(coordinates),direction)))
                        #first resolve overlaps with current gene dict
                        hmmer_hits = resolve_conflicts(hit_dict,minDomSize=minDomSize,verbose = verbose)
                        # check if any hits pass the filter
                        if len(hmmer_hits) > 0:
                            for hmmerHit in hmmer_hits.values():
                                protein.add_hit(anotID,hmmerHit[0],hmmerHit[1])
                        hit_dict = {}
                    current_gene_id = gene_id
                    hmmer_hit = hit[3]
                    hmm_coverage_start = int(hit[15])
                    hmm_coverage_end = int(hit[16])
                    hmm_length = int(hit[2])
                    hmm_coverage = (hmm_coverage_end-hmm_coverage_start+1.)/hmm_length
                    gene_hmm_start = int(hit[17])
                    gene_hmm_end = int(hit[18])
                    hit_dict[(gene_hmm_start,gene_hmm_end)] = (hmmer_hit,score,hmm_coverage)
    return proteinDict
def export_hit_summary(multiClusterDict,outfile,hitDictID,maxJump = 100,minClusterSize = 2, hitsToConsider = set(),
                       hitsToIgnore = set(),writeFile=False):
    # First Unpack Cluster Analysis to Filter Eligible Clusters
    filtered_clusters = SortedListWithKey(key=lambda x: -clusterHits)
    for (species,species_clusters) in multiClusterDict.items():
        for cluster in species_clusters:
            clusterHits = len(cluster)
            if hitsToConsider:
                clusterHits = sum([1 for protein in cluster if len(hitsToConsider & protein.hit_dict[hitDictID].hits) > 0]) - \
                    sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])
            else:
                clusterHits = len(cluster) - \
                    sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])
            # Check size first
            if clusterHits >= minClusterSize:
                proteinIdx = [protein.idx for protein in cluster]
                proteinIdxDiff = [abs(j-i) for i,j in zip(proteinIdx, proteinIdx[1:])]
                proteinMaxJump = max(proteinIdxDiff)
                proteinMaxIdx = proteinIdxDiff.index(proteinMaxJump)
                # if the max gap happens at the start or end of the cluster see if removing that protein will have
                # it fit the threshhold
                while clusterHits > minClusterSize and (max(proteinIdxDiff) > maxJump) and \
                        (proteinMaxIdx in [0,len(proteinIdxDiff)-1]):

                    if proteinMaxIdx == 0:
                        cluster.pop(proteinMaxIdx)
                    else:
                        cluster.pop(proteinMaxIdx+1)

                    if hitsToConsider:
                        clusterHits = sum([1 for protein in cluster if len(hitsToConsider & protein.hit_dict[hitDictID].hits) > 0]) - \
                                    sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])
                    else:
                        clusterHits = len(cluster) - \
                                            sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])

                    proteinIdx = [protein.idx for protein in cluster]
                    proteinIdxDiff = [abs(j-i) for i,j in zip(proteinIdx, proteinIdx[1:])]
                    proteinMaxJump = max(proteinIdxDiff)
                    proteinMaxIdx = proteinIdxDiff.index(proteinMaxJump)

                if max(proteinIdxDiff) < maxJump:
                    if hitsToConsider:
                        clusterHits = sum([1 for protein in cluster if len(hitsToConsider & protein.hit_dict[hitDictID].hits) > 0]) - \
                                    sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])
                    else:
                        clusterHits = len(cluster) - \
                                            sum([1 for protein in cluster if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) > 0])
                    print(cluster)
                    filtered_clusters.add(cluster)
    print("Found %i clusters" % len(filtered_clusters))
    if writeFile:
        with open(outfile,'w') as outHandle:
            outHandle.write('Protein Hit\tProtein Idx\tCluster Hit\n')
            ctr = 1
            for cluster in filtered_clusters:
                outHandle.write('# Species: %s\tNumber Proteins: %i\tDNA Size: %i\n' % (cluster[0].species,len(cluster),cluster.size()))
                outHandle.write('# Protein Hit\tProtein Idx\tCluster Hit\n')
                for protein in cluster:
                    if len(hitsToIgnore & protein.hit_dict[hitDictID].hits) >= 1:
                        outHandle.write('%s\t%i\t##%s\n'% (protein.name,protein.idx,list(protein.hit_dict[hitDictID].hits & hitsToIgnore)[0][0]))
                    elif len(hitsToConsider & protein.hit_dict[hitDictID].hits) >= 1:
                        outHandle.write('%s\t%i\t**%s\n'% (protein.name,protein.idx,list(protein.hit_dict[hitDictID].hits & hitsToConsider)[0][0]))
                    elif len(protein.hit_dict[hitDictID].hits) >= 1:
                            outHandle.write('%s\t%i\t%s\n'% (protein.name,protein.idx,list(protein.hit_dict[hitDictID].hits)[0][0]))
                    else:
                        outHandle.write('%s\t%i\tNo Hits\n'% (protein.name,protein.idx))

    return filtered_clusters
def add_sequences(fastaFile,proteinDict,speciesFilter = None):
    '''
    :param proteinDict: {'protein_ID':protein}
    :param fastaFile: fasta file of sequences with header as follows >speciesID_record|nucl pos|direction|int ID|protein ID
    :return: proteinDict with sequences added to protein object
    '''
    for entry in SeqIO.parse(fastaFile,'fasta'):
        entry_parse = entry.name.split('|')
        species_id = entry_parse[0]
        if ((not speciesFilter) or (species_id in speciesFilter)):
            coordinates = [int(x) for x in entry_parse[1].split('-')]
            direction = entry_parse[2]
            geneIntID = entry_parse[3]
            geneProtID = entry_parse[4]
            geneIntIdx = int(geneIntID.split('_')[-1])
            protein = proteinDict.setdefault(geneIntID,Protein(species_id,geneProtID,geneIntID,geneIntIdx,(tuple(coordinates),direction)))
            protein.sequence = str(entry.seq[:-1])
    return proteinDict

def merge_annotations(protein,mergedID,anot1,anot2,minDomSize,delOrig = False):
    merged = SortedDict()
    # added code for no hits
    if anot1 in protein.annotations.keys():
        merged.update(deepcopy(protein.annotations[anot1].items()))
    if anot2 in protein.annotations.keys():
        merged.update(deepcopy(protein.annotations[anot2].items()))
    merged = resolve_conflicts(merged,minDomSize=minDomSize)
    if len(merged) > 0:
        protein.annotations[mergedID] = merged
    else:
        protein.annotations[mergedID] = SortedDict({(1,protein.size()):('UNK_%.3i' %
                                                                          round(protein.size()/float(minDomSize)),-1,1.)})
    if delOrig:
        if anot1 in protein.annotations.keys():
            del protein.annotations[anot1]
        if anot2 in protein.annotations.keys():
            del protein.annotations[anot2]

def addUnkDomains(protein,anotID,minDomSize):
    # Edge case for start and end windows
    windows = protein.annotations[anotID].keys()
    if calculate_window((1,windows[0][0]-1)) >= minDomSize:
        protein.annotations[anotID][(1,windows[0][0]-1)] = ('UNK_%.3i' % round(calculate_window((1,windows[0][0]-1))/minDomSize),
                                                            1,1.)
    if calculate_window((windows[-1][-1]+1,protein.size())) >= minDomSize:
        protein.annotations[anotID][(windows[-1][-1]+1,protein.size())] = ('UNK_%.3i' % round(calculate_window((windows[-1][-1]+1,protein.size()))/minDomSize),
                                                            1,1.)
    # First Line: set up unknown domain similar to domain format
    # Second Line: to get the difference of the coordinates
    # Third Line: Check if window in between hits meets the minimum domain size
    iterator = (((firstWindow[1]+1,secondWindow[0]-1),('UNK_%.3i' % round(abs(secondWindow[0]-firstWindow[1]-1)/minDomSize),1,1.))
                for firstWindow,secondWindow in zip(windows,windows[1:])
                if abs(secondWindow[0]-firstWindow[1]-1) >= minDomSize)
    protein.annotations[anotID].update(iterator)
    return protein

def addNoHitsCluster(cluster,anotID,minDomSize):
    for protein in cluster:
        if anotID in protein.annotations.keys():
            addUnkDomains(protein,anotID,minDomSize)
        else:
            protein.annotations[anotID] = SortedDict({(1,protein.size()):('UNK_%.3i' %
                                                                          round(protein.size()/float(minDomSize)),-1,1.)})
    return cluster


def retDomStrings(cluster,anot_key,delim,mergeSameDir=True):
    '''
    :param cluster: cluster of proteins
    :param mergeSameDir: will return a chain of annotated IDs
    :return:
    '''
    if not mergeSameDir:
        for protein in cluster:
            yield delim.join(hit[0] for hit in protein.annotations[anot_key].itervalues())
    else:
        directions = [protein.location[1] for protein in cluster]
        # protein.annotations.get(anot_key,{protein.location[0]:('NO_HIT',-1,1)}).itervalues()
        annotations = ((hit[0] for hit in protein.annotations[anot_key].itervalues())
                       for protein in cluster)
        contig_groups = groupby(iter(zip(count(0),annotations)),lambda x: directions[x[0]])
        for direction,contig_group in contig_groups:
            if direction == '+':
                yield delim.join(map(lambda x: delim.join(x[1]),list(contig_group)))
            else:
                yield delim.join(reversed(map(lambda x: delim.join(x[1]),list(contig_group))))

def retDomStringsProt(cluster,anot_key,delim,mergeSameDir=True):
    '''
    :param cluster: cluster of proteins
    :param mergeSameDir: will return a chain of annotated IDs
    :return:
    '''
    if not mergeSameDir:
        for protein in cluster:
            yield delim.join(hit[0] for hit in protein.annotations[anot_key].itervalues())
    else:
        directions = [protein.location[1] for protein in cluster]
        annotations = ((hit[0] for hit in protein.annotations[anot_key].itervalues())
                       for protein in cluster)
        proteinNames = (protein.name for protein in cluster)
        contig_groups = groupby(iter(zip(count(0),annotations,proteinNames)),lambda x: directions[x[0]])
        for direction,contig_group in contig_groups:
            hitProtPairs = map(lambda x: zip(x[1],repeat(x[2])),list(contig_group))
            if direction == '+':
                yield map(lambda x: delim.join(x),zip(*map(lambda y: map(lambda z:delim.join(z),zip(*y)),hitProtPairs)))
            else:
                yield map(lambda x: delim.join(x),zip(*reversed(map(lambda y: map(lambda z:delim.join(z),zip(*y)),hitProtPairs))))


def filter_hmmsearch_hmmStartEnd(path,hmmstartUpperBound,hmmendLowerBound,outfile,hmmcoverageFilter=0,eval_cutoff=1e-5,cutoff_score=25):
    # delete any existing file
    outHandle = open(outfile,'wb')
    with open(path) as hmmfile:
        for line in hmmfile:
            if line[0] == '#':
                pass
            else:
                hit = line.split()
                score = float(hit[13])
                evalue = float(hit[12])
                if score >= cutoff_score  and evalue <= eval_cutoff:
                    #check if you're done parsing the hits of the old gene and annotate the gen
                    hmm_coverage_start = int(hit[15])
                    hmm_coverage_end = int(hit[16])
                    hmm_length = int(hit[5])
                    hmm_coverage = (hmm_coverage_end-hmm_coverage_start+1.)/hmm_length
                    if hmm_coverage_start < hmmstartUpperBound and hmm_coverage_end > hmmendLowerBound and hmm_coverage > hmmcoverageFilter:
                        with open(outfile,'ab') as outHandle:
                            outHandle.write(line)

def predictDisorder(protein,pathToIUPRED):
    '''
    Given a protein with a sequence will use IUPred to get a per residue probability for IUPred to add in annotations
    requires IUPred installer
    :param protein:
    :param pathToIUPRED:
    :return: protein with IUPred annotation added
    '''
    if protein.sequence == '':
        print("Need sequence information for protein to run this")
        return protein
    else:
        currentDir = os.getcwd()
        tmpFastaFile = os.getcwd() + os.sep + protein.name + ".fasta"
        tmpRecord = SeqRecord(Seq(protein.sequence),id=protein.fastaID)
        SeqIO.write(tmpRecord,open(tmpFastaFile,'w'),'fasta')

        os.chdir(pathToIUPRED)
        IUPredout = subprocess.Popen('./iupred %s short' % tmpFastaFile,shell=True,stdout=subprocess.PIPE)
        protein.annotations['IUPRED'] = []
        for line in IUPredout.stdout:
            if line[0] == '#':
                pass
            else:
                protein.annotations['IUPRED'].append(float(line.strip().split()[-1]))

        assert len(protein.annotations['IUPRED']) == len(protein.sequence)

        os.chdir(currentDir)
        os.remove(tmpFastaFile)

        return protein
