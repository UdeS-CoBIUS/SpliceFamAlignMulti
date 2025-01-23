#!/usr/bin/python3
#-*- coding: utf-8 -*-

"""

``utils_reduce.py`` **module description**:

This module contains functions to: 
- create the transcribed sequence of each gene
- compute the MSpA using the transcribed sequences instead of transcript sequences

.. moduleauthor:: Aida Ouangraoua, Chakirou Alabani

2024

"""


import os
import argparse
import glob
import multiprocessing
from functools import partial
from contextlib import contextmanager
from multiprocessing import Pool
import time

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO


# compute the transcribed sequence of each gene
def compute_reducedsource_from_files(sourcedata,targetdata):
    sourcedata_r,geneexon_r,cdsexon_r,cds2geneid_r,cds2geneexon_r,gene2cds_r,geneexon = [],{},{},{},{},{},{}
    geneid2seq = {}
    for i in range(len(targetdata)):
        geneid,geneseq = targetdata[i]
        geneid2seq[geneid] = geneseq
        cdsid = geneid+"_r"
        sourcedata_r.append([cdsid,"",geneid,[]]) 
        geneexon_r[geneid] = []
        geneexon[geneid] = []
        cds2geneid_r[cdsid] = geneid
        gene2cds_r[geneid] = [cdsid]
        cdsexon_r[cdsid] = []
        cds2geneexon_r[cdsid] = {}

    
    for i in range(len(sourcedata)):
        cdsid,seq,cdsgeneid,exonlist = sourcedata[i]
        geneexon[cdsgeneid] += [exon[2:] for exon in exonlist]

    for geneid in geneexon.keys():
        cdsid = geneid+"_r"
        geneexon[geneid].sort()
        for i in range(len(geneexon[geneid])-1):
            if(geneexon[geneid][i+1][0] <= geneexon[geneid][i][1]):
                geneexon[geneid][i+1] = [min(geneexon[geneid][i][0],geneexon[geneid][i+1][0]),max(geneexon[geneid][i][1],geneexon[geneid][i+1][1])]
            else:
                geneexon_r[geneid].append(geneexon[geneid][i])
        geneexon_r[geneid].append(geneexon[geneid][-1])
        
    for i in range(len(sourcedata_r)):
        cdsid,seq,cdsgeneid,exonlist = sourcedata_r[i]
        cstop = 0
        for j in range(len(geneexon_r[cdsgeneid])):
            gstart,gstop = geneexon_r[cdsgeneid][j]
            seq +=  geneid2seq[cdsgeneid][gstart:gstop]
            cstart = cstop
            cstop = cstart +(gstop-gstart)
            exonlist.append([cstart,cstop,gstart,gstop])
            cdsexon_r[cdsid].append([cstart,cstop])
            cds2geneexon_r[cdsid][cstart] = [gstart,gstop]
        sourcedata_r[i][1] = seq
        sourcedata_r[i][3] = exonlist
    
    return sourcedata_r,geneexon_r,cdsexon_r,cds2geneid_r,cds2geneexon_r,gene2cds_r

# compute new results (segment matches) from initial pairwise comparison    
def comparisonresults2comparisonresults_idty(comparisonresults,extendedsourcedata,targetdata,nbinitialsource):
    
    #file format
    #len(comparisonresults)=len(comparisonresults_idty)=nb_genes;
    #len(comparisonresults[i])=len(comparisonresults_idty[i])=nb_pairwise_aln_with_genei;
    #comparisonresults[i][j]=[0,aln,0,0,0]; jth alignment with genei
    #aln =[block1,block2,...]; blocki =[cdsstart,cdsend, genestart,geneend]
    
    #comparisonresults_idty[i][j]= [pid_block1,pid_block2,...]; pid of blocks

    comparisonresults_new = []
    comparisonresults_idty = []
    
    datalist = []
    cdsexon = {}
    j= 0
    for gene in targetdata:
        geneid,geneseq = gene
        
        for i in range(nbinitialsource):
            cds = extendedsourcedata[i]
            cdsid,cdsseq,cdsgeneid ,null = cds
            
            status, blocklist,splicing_sites,ttargetcds, texisting = comparisonresults[j][i]
            datalist.append([cdsid,geneid,cdsseq, cdsgeneid, geneseq, blocklist, status])
            if(geneid==cdsgeneid):
                cdsexon[cdsid]=blocklist
        j += 1
    for i in range(len(datalist)):
        datalist[i].append(cdsexon[datalist[i][0]])
            

    p = Pool(multiprocessing.cpu_count())
    pool_results = p.map(partial(pool_compute_segments), datalist)

        
    results = {}
    for item in pool_results:
        results[item[0] + item[1]] = item[2]
            
    #for data in datalist:
        #item = pool_compute_idty(data)
        #results[item[0] + item[1]] = item[2]
        
    for gene in targetdata:
        geneid,geneseq = gene
        comparisonresults_new.append([])
        comparisonresults_idty.append([])
        for i in range(nbinitialsource):
            comparisonresults_new[-1].append([0,[],0,0,0])
            comparisonresults_idty[-1].append([])
            cds = extendedsourcedata[i]
            cdsid,cdsseq,cdsgeneid ,null = cds
            for match in results[cdsid + geneid]:
                
                genestart, cdsstart, length, pid = match
                comparisonresults_new[-1][i][1].append([cdsstart,cdsstart+length,genestart,genestart+length])
                comparisonresults_idty[-1][i].append(float(pid))

    return comparisonresults_new, comparisonresults_idty
            

def pool_compute_segments(data):
    cdsid,geneid,cdsseq, cdsgeneid,geneseq, blocklist, status,cdsexon = data
    results = []
    segment_matches = []
    
    cdslength = len(cdsseq)
    segments = compute_segments(cdsid,cdsgeneid,geneid,cdsseq, geneseq, blocklist,cdsexon)
    segment_matches += segments

    return [cdsid,geneid,segment_matches]



def compute_segments(cdsid, cdsgeneid,geneid, cds, gene, blocklist,cdsexon):

    segment_matches = []
    cds_len = len(cds)

    if(len(blocklist) > 0):
   
        for i in range(len(blocklist)):
            block = blocklist[i]
            segments = compute_aln(cdsid, cdsgeneid, geneid, cds, gene, block,cdsexon)
            segment_matches += segments
        
    return segment_matches

def compute_aln(cdsid, cdsgeneid, geneid, cds, gene, block, cdsexon):
 
    block_qs = block[0]  # query start
    block_qe = block[1]  # query start
    block_ss = block[2]  # subject start
    block_se = block[3]  # subject end
    # block_identity = "%.2f" % (compute_block_identity(cds, gene, block))
    gene_ = gene[block_ss:block_se]
    cds_ = cds[block_qs:block_qe]

    sequence1 = ""
    sequence2 = ""
    block_identity = 0.0
    if(len(cds_) == len(gene_)):
        sequence1 = gene_
        sequence2 = cds_
    elif(len(cds_) == 0):
        sequence1 = gene_
        sequence2 = '-' * len(sequence1)
    elif(len(gene_) == 0):
        sequence2 = cds_
        sequence1 = '-' * len(sequence2)
    else:
        alignment = pairwise2.align.globalms(gene_, cds_, 2, 0, -10, -1)
        sequence1, sequence2 = alignment[0][0], alignment[0][1]
    aln_length = len(sequence1)

    block_identity = "%.2f" % (1.0 * computeAlignmentPercentIdentity(sequence1, sequence2) / 100)
    
    segment_matches = compute_segment_matches(sequence1, sequence2, block_ss, block_qs, block_identity)

    segment_matches = cut_segment_matches(segment_matches, cdsexon)
    return segment_matches

def compute_segment_matches(sequence1, sequence2, block_ss, block_qs, block_identity):
    segment_matches = []
    i1 = 0
    i2 = 0
    s1 = -1
    s2 = -1
    e1 = -1
    e2 = -1
    for i in range(len(sequence1)):
        if (sequence1[i] != '-' and sequence2[i] != '-'):
            if(s1 == -1):
                s1 = i1
                s2 = i2
            i1 += 1
            i2 += 1
            if(i == len(sequence1) - 1):
                e1 = i1
                e2 = i2
                segment_matches.append([s1+block_ss,s2+block_qs,e1-s1, block_identity])
        else:
            if(s1 != -1):
                e1 = i1
                e2 = i2
                segment_matches.append([s1+block_ss,s2+block_qs,e1-s1, block_identity])
                s1 = -1
                s2 = -1
                e1 = -1
                e2 = -1
            if(sequence1[i] != '-'):
                i1 += 1
            if(sequence2[i] != '-'):
                i2 += 1
    return segment_matches


def cut_segment_matches(segment_matches,cdsexon):
    segment_matches_ = []
    for match in segment_matches:
        genestart,cdsstart,length,pid = match
        cut=[cdsstart]
        for exon in cdsexon:
            if (cdsstart < exon[0] and exon[0] < cdsstart+length):
                cut.append(exon[0])
        cut.append(cdsstart+length)
        if(len(cut)>2):
            print(cut)
        for i in range(1,len(cut)):
            cstart = cut[i-1]
            clength = cut[i]-cut[i-1]
            diffstart=cstart-cdsstart
            segment_matches_.append([genestart+diffstart,cstart,clength,pid])
    return segment_matches_

def computeAlignmentPercentIdentity(sequence1, sequence2):
    """
    This function aim to return the percentage of identity of 2 sequences
    sequence1: string
        sequence 1
    sequence2: string
        sequence 2

     Returns
     -------
     aln_identity: float
        % of identity between the 2 sequences
    """

    aln_identity = 0.0
    match = 0
    length = len(sequence1)

    for i in range(length):

        if (i < len(sequence2) and sequence1[i] == sequence2[i]):
            match += 1
    if length == 0:
        length = 1
    aln_identity = 100.0 * match / length
    return aln_identity

