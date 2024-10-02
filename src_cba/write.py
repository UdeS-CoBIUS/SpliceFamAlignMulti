#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``write.py`` **module description**:

This module is the module that formats and writes the results in output files.

.. moduleauthor:: Aida Ouangraoua

2019

"""

import os
import multiprocessing
from functools import partial
from contextlib import contextmanager
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

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
