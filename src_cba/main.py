#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``main.py`` **module description**:

This module is the main module that 
- parses the input arguments,
- read the input pairwise aligned segment  file
- compute multiple alignment using a consistency-based approach
- computes ortho groups based on multiple alignment
- writes the results in output files

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

from skbio import DistanceMatrix
from skbio.tree import nj

import sys
sys.path.insert(0, os.path.abspath('src_python3/'))
from compute_alignments import *

from compare import *
#from compute_orthology_multi import *


def get_data_from_files(args):
    #Output file format
    #len(sourcedata)=nb_cds; sourcedata[i]=[cdsid,cdsseq,cdsgeneid,exonlist];
    #exonlist =[exon1,exon2,...]; exoni =[cdsstart,cdsend, genestart,geneend]
    #len(targetdata)=nb_genes; targetdata[i]=[geneid,geneseq]
    
    sourcedata = []
    targetdata = []
    sourcefile = args.sourceFile
    if (sourcefile == None):
        print("Argument -sf <sourcefilename> is required")
    targetfile = args.targetFile
    if (targetfile == None):
        print("Argument -tf <targetfilename> is required")
    source2targetfile = args.source2TargetFile
    if (source2targetfile == None):
        print("Argument -s2tf <source2targetfile> is required")

    if(sourcefile != None and targetfile != None and source2targetfile != None):
        for record in SeqIO.parse(sourcefile, "fasta"):
            sourcedata.append([record.id,str(record.seq),"",[]])
        
        for record in SeqIO.parse(targetfile, "fasta"):
            targetdata.append([record.id,str(record.seq)])
        
        source2target = [line.split("\n")[0].split(" ") for line in open(source2targetfile,"r").readlines()]
        for i in range(len(sourcedata)):
            
            sourcedata[i][2] = source2target[i][1]

        sourceexonfile = args.sourceExonFile
        if(sourceexonfile != None):
            file = open(sourceexonfile, "r")
            lines = file.readlines()
            i = 0
            j = 0
            while j < len(lines):
                line = lines[j]
                if(line.startswith('>')):
                    j += 1
                    line = lines[j]
                    exonlist = []
                    while(j < len(lines) and len(lines[j]) > 2 and not(lines[j].startswith('>'))):
                        line = lines[j]
                        tab = line.split("\n")[0].split(" ")
                        exonlist.append([int(x) for x in tab])
                        j += 1
                    sourcedata[i][3] =  exonlist
                    i += 1
                else:
                    j += 1
    return  sourcedata, targetdata


def parseResultFile(resultfile,sourcedata,idty_threshold):
    #Output file format
    #len(comparisonresults)=len(comparisonresults_idty)=nb_genes;
    #len(comparisonresults[i])=len(comparisonresults_idty[i])=nb_pairwise_aln_with_genei;
    #comparisonresults[i][j]=[0,aln,0,0,0]; jth alignment with genei
    #aln =[block1,block2,...]; blocki =[cdsstart,cdsend, genestart,geneend]
    
    #comparisonresults_idty[i][j]= [pid_block1,pid_block2,...]; pid of blocks
    
    #geneexon={}; geneexon[geneid]=[firstexon_cds1,...,lastexon_cds1,firstexon_cds2,...,lastexon_cds2,...]; exon=[start,stop]
    
    #cdsexon={}; cdsexon[cdsid]=[exon1, exon2,...]; exon=[start,stop]

    #cds2gene={}; cds2gene[cdsid] = cdsgeneid

    #cds2geneexon={}; cds2geneexon[cdsid]={}; cds2geneexon[cdsid][cdsexonstart]= correspondinggeneexon

    #gene2cds={}; gene2cds[geneid] = [cdsid1,cdsid2, ....]
 
    file = open(resultfile, "r")
    lines = file.readlines()
    prev_GENE = ""
    liste_GENE    =	[]
    cds2gene = {}
    gene2cds = {}
    cds2geneexon = {}

    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid,null = cds
        cds2gene[cdsid] = cdsgeneid
        
    comparisonresults			=	[]
    comparisonresults_idty		=	[]
    geneexon = {}
    cdsexon = {}
    liste_CDS_GENE		=	[]
    liste_CDS_GENE_idty		=	[]
    i 					=	0
    while i < len(lines):
        line = lines[i]		
        tab = line.split("\t")
        if len(tab) == 5:
            idCDS  		= 	tab[0]
            idGene 		= 	tab[1]
            i += 1
            CDS_GENE = [0, [], 0, 0, 0]
            CDS_GENE_idty = []
            while i< len(lines) and not(lines[i] in ['\n', '\r\n']):
                line = lines[i]	
                tab = line.split("\t")
                if len(tab) >= 8:
                    idCDS			= 	tab[0]
                    idGene			= 	tab[1]
                    length		= 	int(tab[2])
                    beginCDS		= 	int(tab[3])
                    endCDS			= 	int(tab[4])
                    beginGene		=	int(tab[5])
                    endGene			=	int(tab[6])
                    part_CDS_Gene 	= 	[beginCDS, endCDS, beginGene, endGene]
                    part_CDS_Gene_idty = float(tab[7])
                    # if(idGene == cds2gene[idCDS]):
                    #     if(part_CDS_Gene_idty != 1.0):
                    #         print(idCDS,idGene)
                    
                    if(part_CDS_Gene_idty >= idty_threshold and 0 <= beginCDS <= endCDS and 0 <= beginGene <= endGene):
                        CDS_GENE[1].append(part_CDS_Gene)
                        CDS_GENE_idty.append(part_CDS_Gene_idty)
                i += 1

            if(idGene not in geneexon.keys()):
                geneexon[idGene] = []
                gene2cds[idGene]  = []
                
            if(idGene == cds2gene[idCDS]):
                gene2cds[idGene].append(idCDS)
                cdsexon[idCDS] = [x[:2] for x in CDS_GENE[1]]
                geneexon[idGene] += [x[2:] for x in CDS_GENE[1]]
                cds2geneexon[idCDS] = {}
                for x in CDS_GENE[1]:
                    cds2geneexon[idCDS][x[0]] = x[2:]
                
            if(idGene == prev_GENE):
                liste_GENE.append(CDS_GENE)
                liste_GENE_idty.append(CDS_GENE_idty)
            else:
                if(len(liste_GENE) != 0):
                    comparisonresults.append(liste_GENE)
                    comparisonresults_idty.append(liste_GENE_idty)
                liste_GENE    =	[CDS_GENE]
                liste_GENE_idty    =	[CDS_GENE_idty]
            prev_GENE = idGene
        else:
            i += 1

    comparisonresults.append(liste_GENE)
    comparisonresults_idty.append(liste_GENE_idty)
   
    return comparisonresults,comparisonresults_idty,geneexon,cdsexon,cds2gene,cds2geneexon,gene2cds

def compute_tree(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,nbinitialsource):
    rank = {}
    cds2gene = {}
    coverage_matrix = []
    identity_matrix = []
    distance_matrix = []
    list_geneid = []
    for i in range(len(targetdata)):
        coverage_matrix.append([])
        identity_matrix.append([])
        distance_matrix.append([])
        for j in range(len(targetdata)):
            coverage_matrix[i].append([])
            identity_matrix[i].append([])
            distance_matrix[i].append(0)
            
    i=0        
    for gene in targetdata:
        geneid,null = gene
        list_geneid.append(geneid)
        rank[geneid] = i
        i += 1

    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,null,cdsgeneid,null = cds
        cds2gene[cdsid] = cdsgeneid
        j += 1

    i = 0
    for gene in targetdata:
        geneid,null = gene
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            null,blocklist,null,null,null = comparisonresults[i][j]
            coverage = 0
            identity = 0
            k = 0
            for block in blocklist:
                if(block[1] > block[0]):
                    coverage += block[1]-block[0]
                    identity += comparisonresults_idty[i][j][k] * (block[1]-block[0])
                k += 1
            if(coverage != 0):
                identity = 1.0 * identity /coverage
            coverage = 1.0 * coverage /len(cdsseq)
            #print(identity)
            
            coverage_matrix[rank[geneid]][rank[cdsgeneid]].append(coverage)
            coverage_matrix[rank[cdsgeneid]][rank[geneid]].append(coverage)
            identity_matrix[rank[geneid]][rank[cdsgeneid]].append(identity)
            identity_matrix[rank[cdsgeneid]][rank[geneid]].append(identity)
            j += 1
        i += 1
    for i in range(len(targetdata)):
        for j in range(len(targetdata)):
            cover = 0.0
            idty = 0.0
            if(len(coverage_matrix[i][j]) > 0):
               cover = sum(coverage_matrix[i][j])/len(coverage_matrix[i][j])
            if(len(identity_matrix[i][j]) > 0):
               idty = sum(identity_matrix[i][j])/len(identity_matrix[i][j])
            distance_matrix[i][j]= 1.0 - (cover*idty)
            if(i==j and (distance_matrix[i][j] != 0)):
                #print(targetdata[i][0])
                distance_matrix[i][j] = 0
    
    #print(distance_matrix)
    if(len(targetdata) >= 3):
        distance_matrix = DistanceMatrix(distance_matrix, list_geneid)
        tree = nj(distance_matrix)
        tree = tree.root_at_midpoint()
        tree = str(tree)
        tree = tree.split("root")[0]+";"
    elif(len(targetdata) == 2):
        tree = "("+list_geneid[0]+","+list_geneid[1]+");"
    else:
        tree = list_geneid[0]+";"
    
    return tree

#############################
##  OUTPUT FILE WRITTING #####
############################

def write_output_files(extendedsourcedata,targetdata,nbinitialsource,mblocklist,outputprefix):
    
    macroalignment_buffer = ""

    i = 0
    
    for mblock in mblocklist:
        macroalignment_buffer += ">"+ str(i)+"\n"
        for gene in targetdata:
            geneid,geneseq = gene
            if geneid in mblock.keys():
                macroalignment_buffer += geneid + ":" + str(mblock[geneid][0])+"-"+str(mblock[geneid][1])+"\n"
            else:
                macroalignment_buffer +=geneid + ":0-0\n"
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            if cdsid in mblock.keys():
                macroalignment_buffer += cdsid + ":" + str(mblock[cdsid][0])+"-"+str(mblock[cdsid][1])+"\n"
            else:
                macroalignment_buffer += cdsid + ":0-0\n"
        macroalignment_buffer += "\n"
        i+=1
        
    macroalignmentfile = open(outputprefix+"macroalignment.txt","w")
    macroalignmentfile.write(macroalignment_buffer)
    macroalignmentfile.close()


    aln = {}
    all_ids = []
    for gene in targetdata:
        geneid,geneseq = gene
        aln[geneid] = ""
        all_ids.append(geneid)
    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        aln[cdsid] = ""
        all_ids.append(cdsid)

    mblocklistnum = []
    for i in range(len(mblocklist)):
        mblocklistnum.append([i, mblocklist[i]])
    p = Pool(multiprocessing.cpu_count())
    results = p.map(partial(pool_write_microalignment,targetdata=targetdata,extendedsourcedata=extendedsourcedata,nbinitialsource=nbinitialsource,all_ids=all_ids,outputprefix=outputprefix), mblocklistnum)

    #results = []
    #for mblocknum in mblocklistnum:
        #results.append(pool_write_microalignment(mblocknum,targetdata,extendedsourcedata,nbinitialsource,all_ids,outputprefix))
        
    for id in all_ids:
        aln[id] = ""
    for i in range(len(mblocklist)):
        for id in all_ids:
            aln[id] += results[i][id]


    microalignment_buffer = ""
    for  id in aln.keys():
        microalignment_buffer += ">"+str(id) + "\n" + str(aln[id])+"\n"
    microalignmentfile = open(outputprefix+"microalignment.fasta","w")
    microalignmentfile.write(microalignment_buffer)
    microalignmentfile.close()


def pool_write_microalignment(mblocknum,targetdata,extendedsourcedata,nbinitialsource,all_ids,outputprefix):
    aln = {}
    i = mblocknum[0]
    mblock = mblocknum[1]

    present_ids = []
    length = 0
    for gene in targetdata:
        geneid,geneseq = gene
        if geneid in mblock.keys() and mblock[geneid][1] > mblock[geneid][0]:
            present_ids.append(geneid)
            aln[geneid] = geneseq[mblock[geneid][0]:mblock[geneid][1]]
            length = len(aln[geneid])

    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        if cdsid in mblock.keys() and mblock[cdsid][1] > mblock[cdsid][0]:
            present_ids.append(cdsid)
            aln[cdsid] = cdsseq[mblock[cdsid][0]:mblock[cdsid][1]]
            length = len(cdsseq[mblock[cdsid][0]:mblock[cdsid][1]])

    for id in all_ids:
        if(id not in present_ids):
            aln[id] = '-'*length

    return aln

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
        sourcedata_r[i][1] = seq
        sourcedata_r[i][3] = exonlist
    
    return sourcedata_r,geneexon_r,cdsexon_r,cds2geneid_r,cds2geneexon_r,gene2cds_r

#def computeResultFiles(sourcedata,targetdata,outputPrefix):
    
    
#####################
### Main ############

def build_arg_parser():
    parser = argparse.ArgumentParser(description="Transcript multiple aligner")
    parser.add_argument('-idty', '--identityThreshold', help="Identity threshold: real between 0.0 and 1.0 (default = 0.3)", type=float, default = 0.3)
    parser.add_argument('-treef', '--treeFile', help="tree file name")
    parser.add_argument('-sf', '--sourceFile', help="Source file name (required)")
    parser.add_argument('-tf', '--targetFile', help="Target file name (required)")
    parser.add_argument('-s2tf', '--source2TargetFile', help="Source to target file name (required)")

    parser.add_argument('-sef', '--sourceExonFile', help="Source exon file name")

    parser.add_argument('-palnf', '--pairwiseAlnFile', help="Pairwise alignment file name (required)")

    parser.add_argument('-op', '--outputPrefix', help="Output prefix (required)")

    parser.add_argument('-pscore', '--pairwiseScore', help="Pairwise Score: unit or idty")

    parser.add_argument('-cpairwise', '--computePairwise', help="Compute Pairwise: yes or no")
    
    return parser

def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    treefile = args.treeFile
    idty_threshold = float(args.identityThreshold)
    outputprefix = args.outputPrefix
    pairwisealnfile = args.pairwiseAlnFile
    computepairwise = args.computePairwise

    if (outputprefix == None):
        print("Argument -op <outputPrefix> is required")
        return
            
    if (pairwisealnfile == None and computepairwise != "yes"):
        print("Argument -palnf <pairwiseAlnFile> is required;\n Otherwise Argument -cPairwise must be 'yes'")
        return

    print("Retrieving input data...")
    sourcedata,targetdata = get_data_from_files(args)
    nbinitialsource = len(sourcedata)
    tree = ""
    if (treefile != None):
        tree = open(treefile,"r").readline().split("\n")[0]
    mblocklist = []

    print("Computing multiple alignement...")
    if(outputprefix != None):
        comparisonresults,comparisonresults_idty,geneexon,cdsexon,cds2geneid,cds2geneexon,gene2cds = [],[],{},{},{},{},{}
        if(pairwisealnfile != None):
            print("Reading pairwise alignment file...")
            comparisonresults,comparisonresults_idty,geneexon,cdsexon,cds2geneid,cds2geneexon,gene2cds = parseResultFile(pairwisealnfile,sourcedata, idty_threshold)
            if(treefile == None):
                tree = compute_tree(sourcedata,targetdata,comparisonresults,comparisonresults_idty,nbinitialsource)

            mblocklist = compute_msa(sourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,cds2geneid,cds2geneexon,gene2cds,args.pairwiseScore,tree,outputprefix)
            
        elif(computepairwise == "yes" and args.source2TargetFile != None and args.sourceExonFile != None):
            sourcedata_r,geneexon_r,cdsexon_r,cds2geneid_r,cds2geneexon_r,gene2cds_r = compute_reducedsource_from_files(sourcedata,targetdata)
            nbinitialsource_r = len(sourcedata_r)
            print(sourcedata_r, geneexon_r,cdsexon_r,cds2geneid_r,cds2geneexon_r,gene2cds_r)
            comparisonResults_r = spliceAlignment(sourcedata_r, targetdata, 2, 'Yes', 'splign', 'sfa', outputprefix)
            print(comparisonResults)

            #comparisonResults = spliceAlignment(sourceData, targetData, step, compareExon, choice, pairwiseMethod, outputPrefix)
            #comparisonresults_,comparisonresults_idty_r = computeResultFiles(sourcedata_r,targetdata,outputPrefix)
            #if(treefile == None):
            #    tree = compute_tree(sourcedata_r,targetdata,comparisonresults_,comparisonresults_idty_r,nbinitialsource_r)

            mblocklist = compute_msa(sourcedata_r,targetdata,comparisonresults_r,comparisonresults_idty_r,geneexon_r,cdsexon_r,nbinitialsource_r,cds2geneid_r,cds2geneexon_r,gene2cds_r,args.pairwiseScore,tree,outputprefix)

        
        
    write_output_files(sourcedata,targetdata, nbinitialsource,mblocklist,outputprefix)

if __name__ == '__main__':
    main()
