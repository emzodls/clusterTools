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

import clusterAnalysis,clusterDist,clusterGraphics,clusterGenbank
from glob import glob

# Load the results of a phmmer search of refseq bacterial database for cpkG
# Command Used in phmmer search was:
# phmmer --domtblout cpkG_refseq.out --noali -E 1e-5 cpkG.fa /Volumes/Data/refseq/bact_refseq.fa
prots = dict()
prots = clusterAnalysis.parse_hmmsearch_domtbl_anot('cpkG_refseq.out',35,'cpkG',prots)

# This creates a filter for the parser to save memory as it parses through the NRPS annotations
speciesSet = set(x for x.species in prots)
# NRPS_refseq.out is a preprocessed hmmsearch of the refseq database using the NRPS.hmm from antismash
# the speciesFilter argument ensures that the only pieces of DNA that are considered have cpkG homologs
prots = clusterAnalysis.parse_hmmsearch_domtbl_anot('NRPS_refseq.out',35,'NRPS',prots,speciesFilter=speciesSet)
# Proteins from the set of proteins added that have cpkG homology
cpkGProts = set(protein for protein in prots.itervalues() if len(set(['cpkG']) &
                protein.getAnnotations('cpkG')) >=1)
# Proteins from the set that contain both a Thioester reductase and Ketosynthase domain from the hmmsearch
KS_TRprots = set(protein for protein in prots.itervalues()
                 if len(set(['R_DOMAIN','KS_DOMAIN']) & protein.getAnnotations('NRPS')) >=2)
# groups proteins that are within 100kb of each other into cluster objects
putativeClusters = clusterAnalysis.cluster_proteins(prots.itervalues(),100000)
# ignores "clusters" that only have one protein
putativeClusters = {species:[cluster for cluster in clusters if len(cluster) >= 2]
                    for species,clusters in putativeClusters.iteritems()}
putativeClusters = {species:clusters for species,clusters in putativeClusters.iteritems() if clusters}

# creates a python set object for each of the clusters for faster comparison
putativeClusterSets = {species:[(set(protein for protein in cluster),cluster) for cluster in clusters]
                       for species,clusters in putativeClusters.iteritems()}

# Find clusters within a 100kb window that have a cpkG hit and a gene with a KS and TR domain
filteredClusters = dict()
for species,clusters in putativeClusterSets.iteritems():
    for hitProts,cluster in clusters:
        if len(KS_TRprots & hitProts) >= 1 and len(cpkGprots & hitProts) >= 1:
            clusterList = filteredClusters.get(species,[])
            clusterList.append(cluster)
            filteredClusters[species] = clusterList

# Now fetch of the putative clusters for further analysis, window from midpoint of 2 genes is 100kb
clusterGenbank.fetchClusterGenbanks(filteredClusters,100000,writeSummary = True)

# Now you can analyze the putative clusters by feeding the fetched genbanks into antismash and
# analyzing the results once you select the clusters from genbank that you want to further analyze you
# can use clusterTools to do further analysis

clustersToAnalyze = glob('clusters/*.gbk')
clusterGenbank.batch_process(clustersToAnalyze, outputPath = '.', inclNtCDS = False, inclProm = False,
                             singleFile=True)

# Run BLAST comparisons to mibig proteins
'''
makeblastdb -dbtype prot -in allOutput_CDS_prot.fasta
blastp -db allOutput_CDS_prot.fasta -query allOutput_CDS_prot.fasta -outfmt 6 -evalue 1E-5 -max_hsps 1 -out self.tsv
blastp -db ../../../MiBIG/mibig_CDS_prot.fasta -query allOutput_CDS_prot.fasta
-outfmt 6 -evalue 1E-5 -max_hsps 1 -out putative-mibigCmp.tsv -num_threads 4
'''

# Loads the MiBIG biosynthetic gene clusters as clusterTool objects
mibigProts = dict()
mibigProts = clusterAnalysis.add_sequences('clusterCmp.fasta',proteinDict=mibigProts)
print len(mibigProts)
for x in mibigProts.itervalues():
    x.species = x.species.split('.')[0]

mibigClusters = clusterAnalysis.cluster_proteins(mibigProts.itervalues(),1e7)
mibigClusters = {k:v[0] for k,v in mibigClusters.iteritems()}

idxFile = open('../id_part_prod_class_spec.csv')
parse = [line.split(',') for line in idxFile][1:]
mibigIDs = [x[0] for x in parse]
types = [set(map(lambda y: y.strip(),x[3].split('/'))) for x in parse]
products = [x[2] for x in parse]
mibigTypes = dict(zip(mibigIDs,types))
mibigID2Prod = dict(zip(mibigIDs,products))

# Dictionary to store the distances between the putative clusters and each of the mibig clusters
bitDist = dict()

# Load the putative clusters for comparison
prots = dict()
prots = clusterAnalysis.add_sequences('allOutput_CDS_prot.fasta',proteinDict=prots)
prots = clusterAnalysis.parse_allVall_blast_file('self.tsv',evalCutoff=1e-5,proteinDict=prots)
prots = clusterAnalysis.parse_allVall_blast_file('putative-mibigCmp.tsv',evalCutoff=1e-5,proteinDict=prots)

#This assumes there's only one putative cluster per species, you can change this as appropriate
putativeClusters = clusterAnalysis.cluster_proteins(prots.itervalues(),1e7)
putativeClusters = {k:v[0] for k,v in putativeClusters.iteritems()}

# calculate the distances between the putative clusters and the mibig clusters
for cluster in putativeClusters.keys():
    similarClusBit = []
    for key in mibigIDs:
        distBit,pairs = clusterDist.calculateDistBitScore(putativeClusters[cluster],mibigClusters[key],'blast')
        if distBit < 1: similarClusBit.append((distBit,key, mibigID2Prod[key]))
    similarClusBit.sort()
    print "Top 5 Similar Bit Distance Clusters for %s" % cluster
    print similarClusBit[:min(4,len(similarClusBit))]
    bitDist[cluster] = similarClusBit

# Generate the cluster comparison graphics
for clusterID in putativeClusters.keys():
    putativeCluster = putativeClusters[clusterID]
    mibigClusterTag = bitDist[clusterID][0]
    dist,mibigID,product = mibigClusterTag
    mibigCluster = mibigClusters[mibigID]
    dist,pairs = clusterDist.calculateDistBitScore(putativeCluster,mibigCluster,'blast')
    clusterGraphics.generateClusterCompGraphic(putativeCluster,mibigCluster,pairs,'%s_%s_%.2f' % (clusterID,product,dist))