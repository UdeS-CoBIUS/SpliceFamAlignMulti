#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compare.py`` **module description**:

This module is the module that computes a multiple spliced alignment from the set of all pairs of source CDS and target gene comparison. The multiple alignment is a progressive alignment following a guide tree. 

.. moduleauthor:: Abigail Djossou and Aida Ouangraoua

2019

"""

import newick
import numpy
from collections import Counter
from Bio import pairwise2

MIN_FUSION = 0.3
MIN_ADD = 0.3
MIN_IDENTITY1 = 0.5
MIN_IDENTITY2 = 0.3
MIN_DIFFLENGTH = 31

def compute_msa(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,tree,compareExon):
    rank = {} #dict: sequence_id -> index of sequences in matrix comparisonresults (lines:genes; columns:cds)
    list_cds = {} #dict: gene_id -> list of its cdsid
    i = 0
    for gene in targetdata:
        geneid,null = gene
        rank[geneid] = i
        list_cds[geneid] = []
        i += 1
    allcdsseq = {}
    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        rank[cdsid] = j
        allcdsseq[cdsid] = cdsseq
        list_cds[cdsgeneid].append(cdsid)
        j += 1

    #print(tree)
    t = newick.loads(tree)


    mblocklist = compute_msa_recursif(extendedsourcedata,nbinitialsource,t[0],rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon)[0]
    
    return mblocklist 

# recursively compute msa at each node of a tree, from leaves to root
def compute_msa_recursif(extendedsourcedata,nbinitialsource,node,rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon):
    mblocklist = []
    geneidList = []
    #if node is a leaf
    if(node.name != None):
        geneid  = node.name
        i = rank[geneid]
        geneidList = [geneid]
        #compute msa between genes and its cds
        for cdsid in list_cds[geneid]:
            j = rank[cdsid]
            null,blocklist,null,null,null = comparisonresults[i][j]
            for block in blocklist:
                gstart,gend = block[2:]
                k = 0
                while k < len(mblocklist) and mblocklist[k][geneid][1] < gstart:
                    k += 1

                # if new block is after all current mblocks
                if(k == len(mblocklist)):
                    mblocklist.append({})
                    mblocklist[-1][geneid] = block[2:]
                    mblocklist[-1][cdsid] = block[:2]
                # if new block is between mblocks k and k+1 
                elif(gend <= mblocklist[k][geneid][0]):
                    mblocklist.insert(k,{})
                    mblocklist[k][geneid] = block[2:]
                    mblocklist[k][cdsid] = block[:2]
                # if new block == mblock k 
                elif(mblocklist[k][geneid] == block[2:]):
                    mblocklist[k][cdsid] = block[:2]
                # if new block overlaps mblock k 
                else:
                    gstart,gend = mblocklist[k][geneid]
                    mblocklist[k][geneid] = [min(gstart,block[2]),max(gend,block[3])]
                    if(cdsid in mblocklist[k].keys()):
                        cstart,cend = mblocklist[k][cdsid]
                        mblocklist[k][cdsid] = [min(cstart,block[0]),max(cend,block[1])]
                    else:
                        mblocklist[k][cdsid] = block[:2]
        ## merge overlapping mblocks                
        for k in range(len(mblocklist)-1,0,-1):
            if(mblocklist[k-1][geneid][0] <  mblocklist[k][geneid][1] and mblocklist[k][geneid][0] <  mblocklist[k-1][geneid][1]):
                for id in list(set(mblocklist[k-1].keys())&set(mblocklist[k].keys())):
                    mblocklist[k-1][id] = [min(mblocklist[k-1][id][0],mblocklist[k][id][0]),max(mblocklist[k-1][id][1],mblocklist[k][id][1])]
                for id in list(set(mblocklist[k].keys())-set(mblocklist[k-1].keys())):
                    mblocklist[k-1][id] = mblocklist[k][id]
                del(mblocklist[k])

                    
    #if node is an internal
    else:
        # compute msa at left child
        mblocklistLeft,geneidLeft  = compute_msa_recursif(extendedsourcedata,nbinitialsource,node.descendants[0],rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon)
        # compute msa at right child
        mblocklistRight,geneidRight = compute_msa_recursif(extendedsourcedata,nbinitialsource,node.descendants[1],rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon)
        # merge left and right into a single msa
        mblocklist = merge(mblocklistLeft,mblocklistRight,geneidLeft,geneidRight,rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon)
        geneidList = geneidLeft + geneidRight

    #print(node)
    #print(mblocklist)
    #print("check")
    check_order(extendedsourcedata,nbinitialsource,mblocklist)
    check(extendedsourcedata,nbinitialsource,mblocklist)

    return mblocklist, geneidList

# merge two msa on two distinct set of sequences into a single msa
def merge(mblocklistLeft,mblocklistRight,geneidLeft,geneidRight,rank,list_cds,comparisonresults,comparisonresults_idty,geneexon,cdsexon,allcdsseq,compareExon):
    mblocklist = []
    cdsidLeft = []
    cdsidRight = []
    #compute left cdsid list
    for geneid in geneidLeft:
        for cdsid in list_cds[geneid]:
            cdsidLeft.append(cdsid)
    #compute right cdsid list
    for geneid in geneidRight:
        for cdsid in list_cds[geneid]:
            cdsidRight.append(cdsid)

    nbgeneRight = []
    nbgeneLeft = []
    nbcdsRight = []
    nbcdsLeft = []
    addsequenceLeft = []
    addsequenceRight = []
    sequenceLeft = []
    sequenceRight = []
    
    for i in range(len(mblocklistLeft)):
        nbgeneLeft.append(len(list(set(list(mblocklistLeft[i].keys()))&set(geneidLeft))))
        nbcdsLeft.append(len(list(set(list(mblocklistLeft[i].keys()))&set(cdsidLeft))))
        addsequenceLeft.append([0])
        sequenceLeft.append([{}])
    for i in range(len(mblocklistRight)):
        nbgeneRight.append(len(list(set(list(mblocklistRight[i].keys()))&set(geneidRight))))
        nbcdsRight.append(len(list(set(list(mblocklistRight[i].keys()))&set(cdsidRight))))
        addsequenceRight.append([0])
        sequenceRight.append([{}])
    
    fusion_count  = []
    fusion_matrix  = []
    for i in range(len(mblocklistLeft)):
        fusion_matrix.append([])
        fusion_count.append([])
        for j in range(len(mblocklistRight)):
            #nbgenej = len(list(set(list(mblocklistRight[j].keys()))&set(geneidRight)))
            #nbcdsj = len(list(set(list(mblocklistRight[j].keys()))&set(cdsidRight)))
            fusion_matrix[i].append(0)
            fusion_count[i].append([0,0])

    for geneid in geneidLeft:
        i = rank[geneid]
        for cdsid in cdsidRight:
            j = rank[cdsid]
            null,blocklist,null,null,null = comparisonresults[i][j]
            r = 0
            for block in blocklist:
                identity = comparisonresults_idty[i][j][r]
                r += 1
                gstart,gend = block[2:]
                cstart,cend = block[:2]
                mblockOverlapGene = []
                kGene = 0
                while kGene < len(mblocklistLeft):
                    if(geneid in mblocklistLeft[kGene].keys() and (mblocklistLeft[kGene][geneid][1] > gstart and gend > mblocklistLeft[kGene][geneid][0])):
                        overlap = 100.0 * min(mblocklistLeft[kGene][geneid][1] - gstart, gend - mblocklistLeft[kGene][geneid][0])/(mblocklistLeft[kGene][geneid][1]-mblocklistLeft[kGene][geneid][0])
                        mblockOverlapGene.append([kGene,overlap])
                    kGene += 1
                mblockOverlapCds = []
                kCds = 0
                while kCds < len(mblocklistRight):
                    if(cdsid in mblocklistRight[kCds].keys() and (mblocklistRight[kCds][cdsid][1] > cstart and cend > mblocklistRight[kCds][cdsid][0])):
                        overlap = 100.0 * min(mblocklistRight[kCds][cdsid][1] - cstart, cend - mblocklistRight[kCds][cdsid][0])/(mblocklistRight[kCds][cdsid][1]-mblocklistRight[kCds][cdsid][0])
                        mblockOverlapCds.append([kCds,overlap])
                    kCds += 1
                        
                if(len(mblockOverlapGene) == 0):
                    for x in mblockOverlapCds:
                        seq = {}
                        seq[geneid] = block[2:]
                        p = searchSequence(seq, sequenceRight[x[0]])
                        if(p==-1):
                            addsequenceRight[x[0]].append(0)
                            sequenceRight[x[0]].append(seq)
                            p = len(sequenceRight[x[0]])-1
                        addsequenceRight[x[0]][p] += 1
                else:
                    for x in mblockOverlapGene:
                        for y in mblockOverlapCds:
                            fusion_matrix[x[0]][y[0]] += 1
                            fusion_count[x[0]][y[0]][0] += 1

    for geneid in geneidRight:
        i = rank[geneid]
        for cdsid in cdsidLeft:
            j = rank[cdsid]
            null,blocklist,null,null,null = comparisonresults[i][j]
            r = 0
            for block in blocklist:
                identity = comparisonresults_idty[i][j][r]
                r += 1
                gstart,gend = block[2:]
                cstart,cend = block[:2]
                mblockOverlapGene = []
                kGene = 0
                while kGene < len(mblocklistRight):
                    if(geneid in mblocklistRight[kGene].keys() and (mblocklistRight[kGene][geneid][1] > gstart and gend > mblocklistRight[kGene][geneid][0])):
                        overlap = 100.0 * min(mblocklistRight[kGene][geneid][1] - gstart, gend - mblocklistRight[kGene][geneid][0])/(mblocklistRight[kGene][geneid][1]-mblocklistRight[kGene][geneid][0])
                        mblockOverlapGene.append([kGene,overlap])
                    kGene += 1
                mblockOverlapCds = []
                kCds = 0
                while kCds < len(mblocklistLeft):
                    if(cdsid in mblocklistLeft[kCds].keys() and (mblocklistLeft[kCds][cdsid][1] > cstart and cend > mblocklistLeft[kCds][cdsid][0])):
                        overlap = 100.0 * min(mblocklistLeft[kCds][cdsid][1] - cstart, cend - mblocklistLeft[kCds][cdsid][0])/(mblocklistLeft[kCds][cdsid][1]-mblocklistLeft[kCds][cdsid][0])
                        mblockOverlapCds.append([kCds,overlap])
                    kCds += 1


                if(len(mblockOverlapGene) == 0):
                    for y in mblockOverlapCds:
                        seq = {}
                        seq[geneid] = block[2:]
                        p = searchSequence(seq, sequenceLeft[y[0]])
                        if(p==-1):
                            addsequenceLeft[y[0]].append(0)
                            sequenceLeft[y[0]].append(seq)
                            p = len(sequenceLeft[y[0]])-1
                        addsequenceLeft[y[0]][p] += 1

                else:
                    for y in mblockOverlapCds:
                        for x in mblockOverlapGene:
                            fusion_matrix[y[0]][x[0]] += 1
                            fusion_count[y[0]][x[0]][1] += 1


    for i in range(len(mblocklistLeft)):
        for j in range(len(mblocklistRight)):
            if(fusion_matrix[i][j]>0):
                nbpair = nbcdsLeft[i]*nbgeneRight[j]+nbgeneLeft[i]*nbcdsRight[j]
                fusion_matrix[i][j]= 1.0 * fusion_matrix[i][j]/nbpair

    for i in range(len(mblocklistLeft)):
        for j in range(len(addsequenceLeft[i])):
            addsequenceLeft[i][j] = 1.0 * addsequenceLeft[i][j]/nbcdsLeft[i]

    for i in range(len(mblocklistRight)):
        for j in range(len(addsequenceRight[i])):
            addsequenceRight[i][j] = 1.0 * addsequenceRight[i][j]/nbcdsRight[i]

    maxAddLeft = []
    for i in range(len(mblocklistLeft)):
        maxAddLeft.append(numpy.argmax(addsequenceLeft[i]))
    maxAddRight = []
    fusionPairs = []
    
    for i in range(len(mblocklistRight)):
        maxAddRight.append(numpy.argmax(addsequenceRight[i]))

    if(len(mblocklistLeft)>0 and len(mblocklistRight)>0):
        maxFusionLeft = numpy.argmax(fusion_matrix, axis=1)
        maxFusionRight = numpy.argmax(fusion_matrix, axis=0)

        for i in range(len(mblocklistLeft)):
            maxFusioni = fusion_matrix[i][maxFusionLeft[i]]
            for j in range(len(mblocklistRight)):
                maxFusionj = fusion_matrix[maxFusionRight[j]][j]
                if(maxFusioni >= MIN_FUSION and maxFusioni ==  maxFusionj and fusion_matrix[i][j] == maxFusioni and min(fusion_count[i][j]) > 0):
                    fusionPairs.append([i,j,maxFusioni])

                
    nbCompatiblePairs = []
    nbBeforePair = []
    before = 0
    after = 1
    orderLeft = partial_order(mblocklistLeft)
    orderRight = partial_order(mblocklistRight)
    conflictPairs = []
    for k in range(len(fusionPairs)):
        nbCompatiblePairs.append(0)
        nbBeforePair.append([])
    for k in range(len(fusionPairs)):
        ik,jk,maxk = fusionPairs[k]
        for l in range(k+1,len(fusionPairs)):
            il,jl,maxl = fusionPairs[l]
            if((il not in orderLeft[ik][after] or jl not in orderRight[jk][before]) and (il not in orderLeft[ik][before] or jl not in orderRight[jk][after])):
                if (il in orderLeft[ik][after] or jl in orderRight[jk][after]) and all([ik not in orderLeft[x][after] for x in range(il,ik)]) and all([jk not in orderRight[x][after] for x in range(jl,jk)]):
                    nbCompatiblePairs[k] += 1
                    nbCompatiblePairs[l] += 1
                    nbBeforePair[l].append(fusionPairs[k][:2])
                            
                elif(il in orderLeft[ik][before] or jl in orderRight[jk][before])  and all([il not in orderLeft[x][after] for x in range(ik,il)]) and all([jl not in orderRight[x][after] for x in range(jk,jl)]):
                    nbCompatiblePairs[k] += 1
                    nbCompatiblePairs[l] += 1
                    nbBeforePair[k].append(fusionPairs[l][:2])
                else:
                    nbCompatiblePairs[k] += 1
                    nbCompatiblePairs[l] += 1
            else:
                conflictPairs.append([k,l])
    
    while(len(conflictPairs) > 0):
        #print("Conflict: ",[[fusionPairs[k][:2],fusionPairs[l][:2]] for k,l in conflictPairs ])
        nbConflict = {}
        for [k,l] in conflictPairs:
            nbConflict[k]=0
            nbConflict[l]=0
        for [k,l] in conflictPairs:
            nbConflict[k] += 1
            nbConflict[l] += 1
        maxConflict = 0
        for [k,l] in conflictPairs:
            if(nbConflict[k] > maxConflict):
                maxConflict = nbConflict[k]
            if(nbConflict[l] > maxConflict):
                maxConflict = nbConflict[l]
        pairs = []
        for [k,l] in conflictPairs:
            if(k not in pairs and nbConflict[k] == maxConflict):
                pairs.append(k)
            if(l not in pairs and nbConflict[l] == maxConflict):
                pairs.append(l)

        minsupp = fusionPairs[pairs[0]][2]
        delk = pairs[0]
        for k in pairs:
            if(fusionPairs[k][2] <  minsupp):
                minsupp = fusionPairs[k][2]
                delk = k
        pair = fusionPairs[delk][:2]
        #print("Delete:",pair )
        del(fusionPairs[delk])
        del(nbCompatiblePairs[delk])
        del(nbBeforePair[delk])
        for x in range(len(fusionPairs)):
            if(pair in nbBeforePair[x]):
                nbBeforePair[x].remove(pair)

        for i in range(len(conflictPairs)-1,-1,-1):
            if(delk in conflictPairs[i]):
                del(conflictPairs[i])
            else:
                if(conflictPairs[i][0] > delk):
                    conflictPairs[i][0] -= 1
                if(conflictPairs[i][1] > delk):
                    conflictPairs[i][1] -= 1
                    
    locateAddLeft = []
    locateAddRight = []
    listAddLeft = []
    listAddRight = []
    
    for i in range(len(mblocklistLeft)):
        beforei = set([-1])
        afteri = set([len(mblocklistRight)])
        listAddLeft.append([])
        for k in range(len(addsequenceLeft[i])):
            if(addsequenceLeft[i][k] >= MIN_ADD):
                beforek,afterk = locate(sequenceLeft[i][k],mblocklistRight)
                beforei = beforei | beforek
                afteri = afteri | afterk
                listAddLeft[i].append(k)
        locateAddLeft.append([beforei,afteri])

    for i in range(len(mblocklistRight)):
        beforei =  set([-1])
        afteri = set([len(mblocklistLeft)])
        listAddRight.append([])
        for k in range(len(addsequenceRight[i])):
            if(addsequenceRight[i][k] >= MIN_ADD):
                beforek,afterk = locate(sequenceRight[i][k],mblocklistLeft)
                beforei = beforei | beforek
                afteri = afteri | afterk
                listAddRight[i].append(k)
        locateAddRight.append([beforei,afteri])
        
    insertLeft = list(set(range(len(mblocklistLeft)))-set([x[0] for x in fusionPairs]))
    insertRight = list(set(range(len(mblocklistRight)))-set([x[1] for x in fusionPairs]))

    nbCompatibleLeft = []
    nbBeforeLeft = []
    for k in range(len(insertLeft)):
        nbBeforeLeft.append([])
        nbCompatibleLeft.append(0)
    nbCompatibleRight = []
    nbBeforeRight = []
    for k in range(len(insertRight)):
        nbBeforeRight.append([])
        nbCompatibleRight.append(0)

    for i in insertLeft:
        if(len(locateAddLeft[i][before] & locateAddLeft[i][after]) != 0):
            addsequenceLeft[i] = []
            sequenceLeft[i] = []
            locateAddLeft[i] = [set(),set()]
    for i in insertRight:
        if(len(locateAddRight[i][before] & locateAddRight[i][after]) != 0):
            addsequenceRight[i] = []
            sequenceRight[i] = []
            locateAddRight[i] = [set(),set()]

    ki = 0
    for i in insertLeft:
        kj = 0
        for j in insertRight:
            if((i not in locateAddRight[j][after] or j not in locateAddLeft[i][after]) and (i not in locateAddRight[j][before] or j not in locateAddLeft[i][before])):
                if (j in locateAddLeft[i][after] or i in locateAddRight[j][before]):
                    nbCompatibleLeft[ki] += 1
                    nbCompatibleRight[kj] += 1
                    nbBeforeRight[kj].append([i,-1])
                            
                elif(j in locateAddLeft[i][before] or i in locateAddRight[j][after]):
                    nbCompatibleLeft[ki] += 1
                    nbCompatibleRight[kj] += 1
                    nbBeforeLeft[ki].append([-1,j])
                else:
                    nbCompatibleLeft[ki] += 1
                    nbCompatibleRight[kj] += 1
            else:
                nbCompatibleLeft[ki] += 1
                nbCompatibleRight[kj] += 1
                if(maxAddLeft[i] <= maxAddRight[j]):
                    addsequenceLeft[i] = []
                    sequenceLeft[i] = []
                    locateAddLeft[i] = [set(),set()]
                
                    for ij in range(len(insertRight)):
                        if([i,-1] in nbBeforeRight[ij] and i not in locateAddRight[insertRight[ij]][before]):
                            nbBeforeRight[ij].remove([i,-1])
                        if([-1,insertRight[ij]] in nbBeforeLeft[ki] and i not in locateAddRight[insertRight[ij]][after]):
                            nbBeforeLeft[ki].remove([-1,insertRight[ij]])

                    for ijk in range(len(fusionPairs)):
                        if([i,-1] in nbBeforePair[ijk] and i not in orderLeft[fusionPairs[ijk][0]][before]):
                            nbBeforePair[ijk].remove([i,-1])
                        if(fusionPairs[ijk][:2] in nbBeforeLeft[ki] and i not in orderLeft[fusionPairs[ijk][0]][after]):
                            nbBeforeLeft[ki].remove(fusionPairs[ijk][:2])
                else:
                    addsequenceRight[j] = []
                    sequenceRight[j] = []
                    locateAddRight[j] = [set(),set()]
                    for ij in range(len(insertLeft)):
                        if([-1,j] in nbBeforeLeft[ij]  and j not in locateAddLeft[insertLeft[ij]][before]):
                            nbBeforeLeft[ij].remove([-1,j])
                        if([insertLeft[ij],-1] in nbBeforeRight[kj] and j not in locateAddLeft[insertLeft[ij]][after]):
                            nbBeforeRight[kj].remove([insertLeft[ij],-1])

                    for ijk in range(len(fusionPairs)):
                        if([-1,j] in nbBeforePair[ijk] and j not in orderRight[fusionPairs[ijk][1]][before]):
                            nbBeforePair[ijk].remove([-1,j])
                        if(fusionPairs[ijk][:2] in nbBeforeRight[kj] and j not in orderRight[fusionPairs[ijk][1]][after]):
                            nbBeforeRight[kj].remove(fusionPairs[ijk][:2])

            kj += 1
        ki += 1

    for l in range(len(insertLeft)):
        i = insertLeft[l]
        for l2 in range(l+1,len(insertLeft)):
            i2 = insertLeft[l2]
            if(i in orderLeft[i2][before]):
                nbCompatibleLeft[l] += 1
                nbCompatibleLeft[l2] += 1
                nbBeforeLeft[l2].append([i,-1])
            elif(i in orderLeft[i2][after]):
                nbCompatibleLeft[l] += 1
                nbCompatibleLeft[l2] += 1
                nbBeforeLeft[l].append([i2,-1])
            else:
                nbCompatibleLeft[l] += 1
                nbCompatibleLeft[l2] += 1

        for k in range(len(fusionPairs)):
            ik,jk,maxk = fusionPairs[k]
            if((i not in orderLeft[ik][after] or jk not in locateAddLeft[i][after]) and (i not in orderLeft[ik][before] or jk not in locateAddLeft[i][before])):
                if (jk in locateAddLeft[i][after] or i in orderLeft[ik][before]):
                    nbCompatibleLeft[l] += 1
                    nbCompatiblePairs[k] += 1
                    nbBeforePair[k].append([i,-1])
                            
                elif(jk in locateAddLeft[i][before] or i in orderLeft[ik][after]):
                    nbCompatibleLeft[l] += 1
                    nbCompatiblePairs[k] += 1
                    nbBeforeLeft[l].append(fusionPairs[k][:2])
                else:
                    nbCompatibleLeft[l] += 1
                    nbCompatiblePairs[k] += 1
            else:
                addsequenceLeft[i] = []
                sequenceLeft[i] = []
                locateAddLeft[i] = [set(),set()]
                nbCompatibleLeft[l] += 1
                nbCompatiblePairs[k] += 1
                
                for ij in range(len(insertRight)):
                    if([i,-1] in nbBeforeRight[ij] and i not in locateAddRight[insertRight[ij]][before]):
                        nbBeforeRight[ij].remove([i,-1])
                    if([-1,insertRight[ij]] in nbBeforeLeft[l] and i not in locateAddRight[insertRight[ij]][after]):
                        nbBeforeLeft[l].remove([-1,insertRight[ij]])

                for ijk in range(len(fusionPairs)):
                    if([i,-1] in nbBeforePair[ijk] and i not in orderLeft[fusionPairs[ijk][0]][before]):
                        nbBeforePair[ijk].remove([i,-1])
                    if(fusionPairs[ijk][:2] in nbBeforeLeft[l] and i not in orderLeft[fusionPairs[ijk][0]][after]):
                        nbBeforeLeft[l].remove(fusionPairs[ijk][:2])
                        
                if(i in orderLeft[ik][before]):
                    nbBeforePair[k].append([i,-1])
                elif(i in orderLeft[ik][after]):
                    nbBeforeLeft[l].append(fusionPairs[k][:2])
                assert((i not in orderLeft[ik][after] or jk not in locateAddLeft[i][after]) and (i not in orderLeft[ik][before] or jk not in locateAddLeft[i][before]))

    
    for l in range(len(insertRight)):
        j = insertRight[l]
        for l2 in range(l+1,len(insertRight)):
            j2 = insertRight[l2]
            if(j in orderRight[j2][before]):
                nbCompatibleRight[l] += 1
                nbCompatibleRight[l2] += 1
                nbBeforeRight[l2].append([-1,j])
            elif(j in orderRight[j2][after]):
                nbCompatibleRight[l] += 1
                nbCompatibleRight[l2] += 1
                nbBeforeRight[l].append([-1,j2])
            else:
                nbCompatibleRight[l] += 1
                nbCompatibleRight[l2] += 1
        for k in range(len(fusionPairs)):
            ik,jk,maxk = fusionPairs[k]
            if((j not in orderRight[jk][after] or ik not in locateAddRight[j][after]) and (j not in orderRight[jk][before] or ik not in locateAddRight[j][before])):
                if (ik in locateAddRight[j][after] or j in orderRight[jk][before]):
                    nbCompatibleRight[l] += 1
                    nbCompatiblePairs[k] += 1
                    nbBeforePair[k].append([-1,j])
                            
                elif(ik in locateAddRight[j][before] or j in orderRight[jk][after]):
                    nbCompatibleRight[l] += 1
                    nbCompatiblePairs[k] += 1
                    nbBeforeRight[l].append(fusionPairs[k][:2])
                else:
                    nbCompatibleRight[l] += 1
                    nbCompatiblePairs[k] += 1
            else:
                addsequenceRight[j] = []
                sequenceRight[j] = []
                locateAddRight[j] = [set(),set()]
                nbCompatibleRight[l] += 1
                nbCompatiblePairs[k] += 1
                for ij in range(len(insertLeft)):
                    if([-1,j] in nbBeforeLeft[ij]  and j not in locateAddLeft[insertLeft[ij]][before]):
                        nbBeforeLeft[ij].remove([-1,j])
                    if([insertLeft[ij],-1] in nbBeforeRight[l] and j not in locateAddLeft[insertLeft[ij]][after]):
                        nbBeforeRight[l].remove([insertLeft[ij],-1])

                for ijk in range(len(fusionPairs)):
                    if([-1,j] in nbBeforePair[ijk] and j not in orderRight[fusionPairs[ijk][1]][before]):
                        nbBeforePair[ijk].remove([-1,j])
                    if(fusionPairs[ijk][:2] in nbBeforeRight[l] and j not in orderRight[fusionPairs[ijk][1]][after]):
                        nbBeforeRight[l].remove(fusionPairs[ijk][:2])

                if(j in orderRight[jk][before]):
                    nbBeforePair[k].append([-1,j])
                elif(j in orderRight[jk][after]):
                    nbBeforeRight[l].append(fusionPairs[k][:2])
                assert((j not in orderRight[jk][after] or ik not in locateAddRight[j][after]) and (j not in orderRight[jk][before] or ik not in locateAddRight[j][before]))

    listmblocks = []
    for k in range(len(fusionPairs)):
        listmblocks.append([nbBeforePair[k], fusionPairs[k][:2],k])
    for k in range(len(insertLeft)):
        listmblocks.append([nbBeforeLeft[k],[insertLeft[k],-1],k])
    for k in range(len(insertRight)):
        listmblocks.append([nbBeforeRight[k],[-1,insertRight[k]],k])

    listmblocks2 = []
    while(len(listmblocks)>0):
        nbPrec =  len(listmblocks)
        next = []
        for x in listmblocks:
            if(len(x[0])==0):
                next.append(x[1])
                listmblocks.remove(x)
        listmblocks2 += next
        for x in next:
            for y in listmblocks:
                if(x in y[0]):
                    y[0].remove(x)
        if(nbPrec ==  len(listmblocks)):
            candidate = []
            for x in listmblocks:
                if(x[1][1]==-1 and all([y[0]== -1 for y in x[0]])):
                    candidate.append([len(x[0]),x[1]])
                elif(x[1][0]==-1 and all([y[1]== -1 for y in x[0]])):
                    candidate.append([len(x[0]),x[1]])
            candidate.sort()
            retained =  candidate[0]
            for x in listmblocks:
                if(x[1] == retained[1]):
                    beforex = x[0]
                    listmblocks.remove(x)
                    listmblocks2.append(retained[1])
                    for y in listmblocks:
                        if(retained[1] in y[0]):
                            y[0].remove(retained[1])
                    if(retained[1][1] == -1):
                        i = retained[1][0]
                        addsequenceLeft[i] = []
                        sequenceLeft[i] = []
                        for y in beforex:
                            addsequenceRight[y[1]] = []
                            sequenceRight[y[1]] = []
                    if(retained[1][0] == -1):
                        j = retained[1][1]
                        addsequenceRight[j] = []
                        sequenceRight[j] = []
                        for y in beforex:
                            addsequenceLeft[y[0]] = []
                            sequenceLeft[y[0]] = []
    for pair in listmblocks2:
        i,j= pair
        if(i != -1 and j != -1):
            mblocklist.append(mblocklistLeft[i])
            for id in list(mblocklistRight[j].keys()):
                if(id in list(mblocklist[-1].keys())):
                    mblocklist[-1][id] = [min(mblocklist[-1][id][0],mblocklistRight[j][id][0]), max(mblocklist[-1][id][1],mblocklistRight[j][id][1])]
                else:
                    mblocklist[-1][id] = mblocklistRight[j][id]
            #print("Fusion",i,j)
        elif(i != -1):
            mblocklist.append(mblocklistLeft[i])
            #print("CopyLeft",i)
            for k in range(len(addsequenceLeft[i])):
                if(addsequenceLeft[i][k] >= MIN_ADD):
                    idk = list(sequenceLeft[i][k].keys())[0]
                    mblocklist[-1][idk] = sequenceLeft[i][k][idk]
                    #print("AddLeft",i,k)
        elif(j != -1):
            mblocklist.append(mblocklistRight[j])
            #print("CopyRight",j)
            for k in range(len(addsequenceRight[j])):
                if(addsequenceRight[j][k] >= MIN_ADD):
                    idk = list(sequenceRight[j][k].keys())[0]
                    mblocklist[-1][idk] = sequenceRight[j][k][idk]
                    #print("AddRight",j,k)

    mblocklist = merge_overlapping(mblocklist)
    if(compareExon == 'Yes'):
        mblocklist = merge_compatible_unordered(mblocklist,allcdsseq)
        
    return mblocklist
                


def order_mblocklist(mblocklist):
    """
    This function orders blocks in the blocklist by increasing
    order of query start location.
    """
    for i in range(len(mblocklist)):
        for j in range(i+1, len(mblocklist)):
            common_keys = list(set(mblocklist[i].keys())&set(mblocklist[j].keys()))
    for i in range(len(mblocklist)):
        min = i
        minprec = -1
        while(min != minprec):
            minprec = min
            for j in range(i,len(mblocklist)):
                common_keys = list(set(mblocklist[j].keys())&set(mblocklist[min].keys()))
                if(len(common_keys) > 0 and all([mblocklist[j][id][1] <=  mblocklist[min][id][0] for id in common_keys])):
                    min = j
                    break
        if(min != i):
            tmp = mblocklist[min]
            mblocklist[min] = mblocklist[i]
            mblocklist[i] = tmp
    return mblocklist

def merge_overlapping(mblocklist):
    mblocklist = order_mblocklist(mblocklist)
    for i in range(len(mblocklist)-1,0,-1):
        for j in range(i-1,-1,-1):
            common_keys = list(set(mblocklist[i].keys())&set(mblocklist[j].keys()))
            if(len(common_keys) > 0):
                overlap = False
                for id in common_keys:
                    if(mblocklist[i][id][0] <  mblocklist[j][id][1] and mblocklist[j][id][0] <  mblocklist[i][id][1]):
                        overlap = True
                        break
                if(overlap):
                    for id in list(mblocklist[i].keys()):
                        if(id in common_keys):
                           mblocklist[j][id] = [min(mblocklist[j][id][0],mblocklist[i][id][0]), max(mblocklist[j][id][1],mblocklist[i][id][1])]
                        else:
                           mblocklist[j][id] = mblocklist[i][id]
                    del(mblocklist[i])
                    break
    mblocklist = order_mblocklist(mblocklist)

    return mblocklist
                    
def partial_order(mblocklist):
    order = []
    before = 0
    after = 1
    for i in range(len(mblocklist)):
        order.append([set([-1]),set([len(mblocklist)])])
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            if(len(set(mblocklist[i].keys()) & set(mblocklist[j].keys())) > 0):
                order[i][after].add(j)
                order[j][before].add(i)
    #for k in range(len(mblocklist)):
    for k in range(len(mblocklist)):
        for i in range(len(mblocklist)):
            for j in (order[i][before]-set([-1,len(mblocklist)])):
                order[i][before] = (order[i][before])|(order[j][before])
            for j in (order[i][after]-set([-1,len(mblocklist)])):
                order[i][after] = (order[i][after])|(order[j][after])
    return order

def merge_compatible_unordered(mblocklist,allcdsseq):
    before = 0
    after = 1
    order = partial_order(mblocklist)
    to_delete = []
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            if(i in to_delete or j in to_delete):
                continue
            if(not(j in (order[i][before])|(order[i][after]))):
                pos_before = (order[i][before])|(order[j][before])
                pos_after = (order[i][after])|(order[j][after])
                
                if (len(pos_before & pos_after)==0):
                    lengthi,seqi = mblocklength(mblocklist[i],allcdsseq)
                    lengthj,seqj = mblocklength(mblocklist[j],allcdsseq)
                    identity = compute_sequence_identity(seqi,seqj)
                    if(lengthi%3 == lengthj%3 and (identity >= MIN_IDENTITY1 or (identity >= MIN_IDENTITY2 and abs(lengthi-lengthj) <= MIN_DIFFLENGTH))):
                        for id in mblocklist[j].keys():
                            mblocklist[i][id]=mblocklist[j][id]
                        order[i][before] = order[i][before]|order[j][before]
                        order[i][after] = order[i][after]|order[j][after]
                        for k in order[i][before]-set([-1,len(mblocklist)]):
                             order[k][after].add(i)
                        for k in order[i][after]-set([-1,len(mblocklist)]):
                             order[k][before].add(i)
                        to_delete.append(j)
                        #print("merged",lengthi,lengthj, lengthi%3 == lengthj%3, identity)
                    #else:
                        #print("not merged", abs(lengthi - lengthj), lengthi%3 == lengthj%3, identity)

    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])

    mblocklist = order_mblocklist(mblocklist)

    return mblocklist

def mblocklength(mblock,allcdsseq):
    lengths = []
    for id in mblock.keys():
        if(id in list(allcdsseq.keys())):
            lengths.append(mblock[id][1]-mblock[id][0])
    length,count = Counter(lengths).most_common(1)[0]
    for id in mblock.keys():
        if(id in list(allcdsseq.keys()) and (mblock[id][1]-mblock[id][0]) == length):
            seq = allcdsseq[id][mblock[id][0]:mblock[id][1]]
            break
    return length,seq

def searchSequence(seq, listSeq):
    pos = -1
    id = list(seq.keys())[0]
    for i in range(len(listSeq)):
        listid = list(listSeq[i].keys())
        if(len(listid)==1 and id  == listid[0] and seq[id][0]<listSeq[i][id][1] and listSeq[i][id][0]<seq[id][1]):
            pos = i
            listSeq[i][id] = [min(seq[id][0],listSeq[i][id][0]),max(seq[id][1],listSeq[i][id][1])]
            break
    return pos
        
def locate(seq,mblocklist):
    before = set([-1])
    after = set([len(mblocklist)])
    id = list(seq.keys())[0]
    for i in range(len(mblocklist)):
        listid = list(mblocklist[i].keys())
        if(id in listid):
            if(mblocklist[i][id][1]<=seq[id][0]):
                before = before | set([i])
            elif(seq[id][1]<=mblocklist[i][id][0]):
                after = after | set([i])
    return before,after
            

def check(extendedsourcedata,nbinitialsource,mblocklist):
    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        prec = [0,0]
        found = False
        for i in range(len(mblocklist)):
            block = mblocklist[i]
            if(cdsid in list(block.keys())):
                found = True
                #print(cdsid,prec,block[cdsid])
                assert(block[cdsid][0] == prec[1])
                prec = block[cdsid]
        if(found):
            assert(prec[1] == len(cdsseq))
        j += 1

def check_order(extendedsourcedata,nbinitialsource,mblocklist):
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            blocki = mblocklist[i]
            blockj = mblocklist[j]
            comparison = {}
            comparison["<"] = 0
            comparison[">"] = 0
            commonids = set(list(blocki.keys())) & set(list(blockj.keys()))
            for id in commonids:
                if(blocki[id][1] <= blockj[id][0]):
                    comparison["<"] += 1
                elif(blockj[id][1] <= blocki[id][0]):
                    comparison["<"] += 1
            assert(comparison["<"] == 0 or comparison[">"] == 0)
            assert(comparison["<"] + comparison[">"] == len(commonids))
    

def compute_sequence_identity(seq1, seq2):
    sequence1 = ""
    sequence2 = ""

    sequence_identity = 0.0

    if(len(seq1)!=0 and len(seq2)!=0):        
        if(len(seq1)==len(seq2)):
            sequence1 = seq1
            sequence2 = seq2
        else:
            # maximise matches and minimize gaps
            alignment = pairwise2.align.globalms(seq1, seq2,2,0,-5,-1)
            sequence1, sequence2 = alignment[0][0],alignment[0][1]
        
        match = 0
        length = len(sequence1)
   
        for i in range(length):
            if(sequence1[i] == sequence2[i]):
                match += 1
        sequence_identity = 1.0 * match /len(sequence2)

    return sequence_identity
