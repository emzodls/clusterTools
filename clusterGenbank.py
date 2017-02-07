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

__author__ = 'emzodls'

from Bio import SeqIO,Seq
from Bio.SeqRecord import SeqRecord
import sys
from Bio import Entrez
import gzip

def fetchClusterGenbanks(clusterDict,window,targetDir ="./",writeSummary = False):
    db = "nuccore"
    Entrez.email = "some_email@somedomain.com"
    retmax = 10**9
    if writeSummary:
        with open('%s/downloadSummary.tsv'% targetDir,'wb') as outfile:
            outfile.write('# Nucleotide Accession\tSpecies\tProtein Hits\tLocation\tWindow Size\n')
    for species,clusters in clusterDict.iteritems():
        species = species.split('.')[0]
        handle = Entrez.esearch( db=db,term=species,retmax=retmax )
        giList = Entrez.read(handle)['IdList']
        search_handle = Entrez.efetch(db=db, id=giList[0], retmode="xml",strand=1,seq_start=1,seq_stop=2)
        name = Entrez.read(search_handle)[0]["GBSeq_definition"]
        search_handle.close()
        sys.stderr.write( "Found %s GI: %s, Species: %s \n" % (len(giList), ", ".join(giList[:10]),name))
        for idx,cluster in enumerate(clusters):
            clusterMidpoint = sum(cluster.location)/2.
            sys.stderr.write( "Fetching cluster %i \n" % (idx + 1) )
            window_start = int(max(1,clusterMidpoint- window/2))
            window_end = int(clusterMidpoint + window/2)
            cluster_filename = "%s/%s-%i-%i.gbk" % (targetDir,species,window_start,window_end)
            try:
                fetch_handle = Entrez.efetch(db=db, id=giList[0], rettype="gbwithparts",strand=1,seq_start=str(window_start),seq_stop=str(window_end))
                with open(cluster_filename,'wb') as genbank_file:
                    genbank_file.write(fetch_handle.read())
                if writeSummary:
                    proteins = ','.join(protein.name for protein in cluster)
                    with open('%s/downloadSummary.tsv'% targetDir,'ab') as outfile:
                        outfile.write('%s\t%s\t%s\t%i-%i\t%i\n' % (species,name,proteins,window_start,window_end,window_end-window_start))
            except IndexError:
                print "Empty GI list for: %s" % species

def fetchGbNucl(targets,window,targetDir ="./",writeSummary = False):
    db = "nuccore"
    Entrez.email = "some_email@somedomain.com"
    retmax = 10**9
    if writeSummary:
        with open('%s/downloadSummary.tsv'% targetDir,'wb') as outfile:
            outfile.write('# Nucleotide Accession\tSpecies\tLocation\tWindow Size\n')
    for species,location in targets:
        species = species.split('.')[0]
        handle = Entrez.esearch( db=db,term=species,retmax=retmax )
        giList = Entrez.read(handle)['IdList']
        search_handle = Entrez.efetch(db=db, id=giList[0], retmode="xml",strand=1,seq_start=1,seq_stop=2)
        name = Entrez.read(search_handle)[0]["GBSeq_definition"]
        search_handle.close()
        sys.stderr.write( "Found %s GI: %s, Species: %s \n" % (len(giList), ", ".join(giList[:10]),name))
        clusterMidpoint = sum(location)/2.
        window_start = int(max(1,clusterMidpoint- window/2))
        window_end = int(clusterMidpoint + window/2)
        cluster_filename = "%s/%s-%i-%i.gbk" % (targetDir,species,window_start,window_end)
        try:
            handle = Entrez.efetch(db=db, id=giList[0], rettype="gbwithparts",strand=1,seq_start=window_start,seq_stop=window_end)
            with open(cluster_filename,'wb') as genbank_file:
                genbank_file.write(handle.read())
            if writeSummary:
                with open('%s/downloadSummary.tsv'% targetDir,'ab') as outfile:
                    outfile.write('%s\t%s\t%i-%i\t%i\n' % (species,name,window_start,window_end,window_end-window_start))
        except IndexError:
            print "Empty GI list for: %s" % species


def batch_process(genbank_file_list, outputPath = '.',speciesOverride = None,
                  offSet = 0, upstream_size = 200, inclNtCDS = True, inclProm = True,singleFile=False,includeORF = False):
    for genbank_file in genbank_file_list:

        print "Parsing: %s" % genbank_file

        if '.gz' in genbank_file:
            genbank_entries = SeqIO.parse(gzip.open(genbank_file),'genbank')
        else:
            genbank_entries = SeqIO.parse(open(genbank_file),"genbank")
        if speciesOverride:
            species_id = speciesOverride
        else:
            species_id = '.'.join(genbank_file.split('/')[-1].split('.')[0:2])
        if inclProm:
            promoters_outfile_name = '%s/%s_promoters.fasta' % (outputPath, species_id)
            outfile = open(promoters_outfile_name,'w')
            outfile.close()
        if inclNtCDS:
            CDS_nt_outfile_name = '%s/%s_CDS_nt.fasta' % (outputPath, species_id)
            outfile = open(CDS_nt_outfile_name,'w')
            outfile.close()

        if not singleFile:
            CDS_prot_outfile_name = '%s/%s_CDS_prot.fasta' % (outputPath, species_id)

        else:
            CDS_prot_outfile_name = '%s/allOutput_CDS_prot.fasta' % outputPath

        cds_ctr = 0
        entry_ctr = 1
        # See if user wants a different name
        for genbank_entry in genbank_entries:
            if speciesOverride:
                species_id = speciesOverride + '_%i' % entry_ctr
            else:
                species_id = genbank_file.split('/')[-1].split('.')[0] + '_%i' % entry_ctr
                if '' != genbank_entry.name:
                    species_id = genbank_entry.name
                elif '' != genbank_entry.id:
                    species_id = genbank_entry.id

            if includeORF:
                CDS_list = (feature for feature in genbank_entry.features if feature.type == 'CDS' or feature.type == 'ORF')
            else:
                CDS_list = (feature for feature in genbank_entry.features if feature.type == 'CDS')

            for CDS in CDS_list:
                cds_ctr += 1
                direction = CDS.location.strand
                # Ensure that you don't get negative values, Biopython parser will not ignore slices that are greater
                # than the entry so you don't need to worry about the other direction
                internal_id = "%s_CDS_%.5i" % (species_id,cds_ctr)
                protein_id = internal_id

                gene_start = max(0,CDS.location.nofuzzy_start)
                gene_end = max(0,CDS.location.nofuzzy_end)
                genbank_seq = CDS.location.extract(genbank_entry)
                nt_seq = genbank_seq.seq

                # Try to find a common name for the promoter, otherwise just use the internal ID
                if 'protein_id' in CDS.qualifiers.keys():
                    protein_id = CDS.qualifiers['protein_id'][0]
                else:
                    for feature in genbank_seq.features:
                        if 'locus_tag' in feature.qualifiers:
                            protein_id = feature.qualifiers['locus_tag'][0]

                if inclProm:
                    promoter_internal = "%s_prom_%.5i" % (species_id,cds_ctr)

                if 'translation' in CDS.qualifiers.keys():
                    prot_seq = Seq.Seq(CDS.qualifiers['translation'][0])
                    if direction == 1:
                        direction_id = '+'
                        if inclProm:
                            promoter_start = max(0,gene_start-upstream_size)
                            promoter_end = max(0,gene_start)
                            assert promoter_start <= promoter_end
                            promoter = genbank_entry.seq[promoter_start:promoter_end]

                    else:
                        direction_id = '-'
                        if inclProm:
                            promoter_start = max(0,gene_end)
                            promoter_end = max(0,gene_end+upstream_size)
                            assert promoter_start <= promoter_end
                            promoter = genbank_entry.seq[promoter_start:promoter_end].reverse_complement()
                else:

                    if direction == 1:
                        direction_id = '+'

                        # for protein sequence if it is at the start of the entry assume that end of sequence is in frame
                        # if it is at the end of the genbank entry assume that the start of the sequence is in frame
                        if gene_start == 0:
                            if len(nt_seq) % 3 == 0:
                                prot_seq = nt_seq.translate()
                            elif len(nt_seq) % 3 == 1:
                                prot_seq = nt_seq[1:].translate()
                            else:
                                prot_seq = nt_seq[2:].translate()
                        else:
                            prot_seq = nt_seq.translate()

                        if inclProm:
                            promoter_start = max(0,gene_start-upstream_size)
                            promoter_end = max(0,gene_start)
                            assert promoter_start <= promoter_end
                            promoter = genbank_entry.seq[promoter_start:promoter_end]

                    if direction == -1:
                        direction_id = '-'

                        nt_seq = genbank_seq.seq

                        if gene_start == 0:
                                prot_seq = nt_seq.translate()
                        else:
                            if len(nt_seq) % 3 == 0:
                                prot_seq = nt_seq.translate()
                            elif len(nt_seq) % 3 == 1:
                                prot_seq = nt_seq[:-1].translate()
                            else:
                                prot_seq = nt_seq[:-2].reverse_complement().translate()
                        if inclProm:
                            promoter_start = max(0,gene_end)
                            promoter_end = max(0,gene_end+upstream_size)
                            assert promoter_start <= promoter_end
                            promoter = genbank_entry.seq[promoter_start:promoter_end].reverse_complement()
                    # Write Nucl file
                if len(nt_seq) > 0 and inclNtCDS:
                    nucl_entry = SeqRecord(nt_seq,id='%s|%i-%i|%s|%s|%s' % (species_id,offSet+gene_start+1,
                                                                            offSet+gene_end,direction_id,internal_id,protein_id),
                                               description = '%s in %s' % (protein_id,species_id))
                    with open(CDS_nt_outfile_name,'ab') as outfile_handle:
                        SeqIO.write(nucl_entry,outfile_handle,'fasta')

                # Write protein file
                if len(prot_seq) > 0:
                    prot_entry = SeqRecord(prot_seq,id='%s|%i-%i|%s|%s|%s' % (species_id,offSet+gene_start+1,
                                                                            offSet+gene_end,direction_id,internal_id,protein_id),
                                               description = '%s in %s' % (protein_id,species_id))
                    with open(CDS_prot_outfile_name,'ab') as outfile_handle:
                        SeqIO.write(prot_entry,outfile_handle,'fasta')

                # Write Promoters File
                if inclProm and len(promoter) > 0:
                    promoter_entry = SeqRecord(promoter,id='%s|%i-%i|%s|%s|%s_upstream' % (species_id,offSet+promoter_start+1,
                                                                    offSet+promoter_end,direction_id,promoter_internal,protein_id),
                        description = '%i nucleotides upstream of %s in %s' % (upstream_size,protein_id,species_id))
                    with open(promoters_outfile_name,'ab') as outfile_handle:
                        SeqIO.write(promoter_entry,outfile_handle,'fasta')
            entry_ctr +=1

def generateSpeciesDict(accList):
    db = "nuccore"
    Entrez.email = "some_email@somedomain.com"
    retmax = 10**9
    speciesDict = dict()
    for species in accList:
        try:
            handle = Entrez.esearch( db=db,term=species,retmax=retmax )
            giList = Entrez.read(handle)['IdList']
            search_handle = Entrez.efetch(db=db, id=giList[0], retmode="xml",strand=1,seq_start=1,seq_stop=2)
            name = Entrez.read(search_handle)[0]["GBSeq_definition"]
            sys.stderr.write( "Found %s GI: %s, Species: %s \n" % (len(giList), ", ".join(giList[:10]),name))
            search_handle.close()
        except:
            name = 'Error Fetching'
        speciesDict[species] = name
    return speciesDict

def removeUndefSeqs(fileName):
    baseName = fileName.split('/')[-1].split('.fasta')[0]
    fastaFile = SeqIO.parse(fileName,'fasta')
    filtered = (rec for rec in fastaFile if any (ch != 'X' for ch in rec.seq))
    outFile = '%s/%s_clean.fasta' % ('/'.join(fileName.split('/')[:-1]),baseName)
    print outFile
    with open(outFile,'ab') as outFileHandle:
        SeqIO.write(filtered,outFileHandle,'fasta')


def batch_process_wrapper(fileName):
    batch_process([fileName],outputPath='/Volumes/Data/new_refseq/',
                  inclNtCDS=False,inclProm=False)
def fetchSeqFromEntry(entry,seqSet,outputFile):
    if entry.id in seqSet:
        with open(outputFile,'ab') as outfile:
            SeqIO.write(entry,outfile,'fasta')


def fetchSeqFromDB(database,seqSet,outputFile):
    for entry in database:
        fetchSeqFromEntry(entry,seqSet,outputFile)

