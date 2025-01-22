#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compare.py`` **module description**:

This module is the module that launches the comparison for all pairs of source CDS and target gene.

.. moduleauthor:: Aida Ouangraoua, Chakirou Alabani

2024

"""

import newick

# Compute the gene location given the cds location
def cds2genelocation(cdsid,segment,cdsexon,cds2geneexon):
    start,end = segment
    gstart,gend = 0,0
    exonstart = -1
    exonend = -1

    for i in range(len(cdsexon[cdsid])):
        exon = cdsexon[cdsid][i]
        if(exon[0] <= start and start < exon[1]):
            gstart = cds2geneexon[cdsid][exon[0]][0] + (start-exon[0])
            exonstart = exon
        if(exon[0] < end and end <= exon[1]):
            gend = cds2geneexon[cdsid][exon[0]][0] + (end-exon[0])
            exonend = exon

    if(exonstart == -1):
        gstart = exonstart
    if(exonend == -1):
        gend = exonend

    return [gstart,gend]
    
# compute segment matches from pairwise comparison results
def comparisonresults2segmentmatches(targetdata, nbinitialsource, extendedsourcedata,comparisonresults,comparisonresults_idty,cdsexon,cds2geneexon,geneid2index):
    segment_matches_matrix=[]
    for i in range(len(targetdata)):
        segment_matches_matrix.append([])
        for j in range(len(targetdata)):
            segment_matches_matrix[i].append([])
    
    i = 0
    # for each block in each alignment
    for gene in targetdata:
        geneid,geneseq = gene
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            
            null,blocklist,null,null,null = comparisonresults[i][j]
            blocklist_idty = comparisonresults_idty[i][j]
            for k in range(len(blocklist)):
                block = blocklist[k]
                gene_segment = block[2:]
                cdssegment = block[:2]
                geneid1 = cdsgeneid
                geneid2 = geneid
                start1 = cds2genelocation(cdsid,cdssegment,cdsexon,cds2geneexon)[0]
                start2 = gene_segment[0]
                length = gene_segment[1]-gene_segment[0]
                idty = blocklist_idty[k]
                if(geneid2index[geneid1] <= geneid2index[geneid2]):
                    segment_matches_matrix[geneid2index[geneid1]][geneid2index[geneid2]].append([start1,start2,length,idty])
                else:
                    segment_matches_matrix[geneid2index[geneid2]][geneid2index[geneid1]].append([start2,start1,length,idty])
        i += 1
                
    for i in range(len(targetdata)):
        for j in range(i,len(targetdata)):
            segment_matches_matrix[i][j].sort()
            for k in range(len(segment_matches_matrix[i][j])-1,0,-1):
                matchk = segment_matches_matrix[i][j][k]
                matchkm1 = segment_matches_matrix[i][j][k-1]
                if( matchk[:3] == matchkm1[:3]):
                    segment_matches_matrix[i][j][k-1][3] = max(segment_matches_matrix[i][j][k-1][3],segment_matches_matrix[i][j][k][3])
                    del(segment_matches_matrix[i][j][k])
        
    return segment_matches_matrix

#Compute refined boundaries given the segments matches
def compute_refine_boundaries(segment_matches_matrix,gene_list):
    genevertices = []#V
    for geneid in gene_list:
        genevertices.append(set())

    for i in range(len(gene_list)):
        for j in range(i,len(gene_list)):
            for match in segment_matches_matrix[i][j]:
                start1,start2,length,idty= match
                genevertices[i].add(start1)
                genevertices[i].add(start1+length)
                genevertices[j].add(start2)
                genevertices[j].add(start2+length)
        

    for i in range(len(gene_list)):
        for j in range(i,len(gene_list)):
            for match in segment_matches_matrix[i][j]:
                start1,start2,length,idty= match
                boundaries = [[i,start1],[i,start1+length],[j,start2],[j,start2+length]]
                for w in boundaries:
                    genevertices = refineboundaries(w,segment_matches_matrix,genevertices)

    for i in range(len(genevertices)):
        genevertices[i] = list(genevertices[i])
        genevertices[i].sort()
                        
    return genevertices

def refineboundaries(w,segment_matches_matrix,genevertices):
    genevertices_ = genevertices
    geneidxw,boundw = w
    overlappingsegments = []
   
    for i in range(geneidxw):
        for k in range(len(segment_matches_matrix[i][geneidxw])):
            match = segment_matches_matrix[i][geneidxw][k]
            start1,start2,length,idty= match
            if(start2 < boundw and boundw <  start2 + length):
                overlappingsegments.append([i,geneidxw,k])
    for i in range(geneidxw,len(genevertices)):
        for k in range(len(segment_matches_matrix[geneidxw][i])):
            match = segment_matches_matrix[geneidxw][i][k]
            start1,start2,length,idty= match
            if(start1 < boundw and boundw <  start1 + length):
                overlappingsegments.append([geneidxw,i,k])
    

    for s in range(len(overlappingsegments)):
        i,j,k = overlappingsegments[s]
        start1,start2,length,idty= segment_matches_matrix[i][j][k]
        geneidxh,boundh = w
        if(geneidxw == i):
            geneidxh = j
            boundh = start2+(boundw-start1)
        else:
            geneidxh = i
            boundh = start1+(boundw-start2)
                
        if(boundh not in genevertices[geneidxh]):
            genevertices_[geneidxh].add(boundh)
            genevertices_ = refineboundaries([geneidxh,boundh],segment_matches_matrix,genevertices_)
    return genevertices_

# cut the segment matches according to the refined boundaries
def refine_segment_matches_matrix(segment_matches_matrix,genevertices):
    segment_matches_matrix_ = []
    for i in range(len(genevertices)):
        segment_matches_matrix_.append([])
        for j in range(len(genevertices)):
            segment_matches_matrix_[i].append([])

    for i in range(len(genevertices)):
        for j in range(i, len(genevertices)):
            for match in segment_matches_matrix[i][j]:
                start1,start2,length,idty= match
                k = 0;
                while(genevertices[i][k] != start1):
                    k += 1
                while(genevertices[i][k] != start1+length):
                    sstart1 = genevertices[i][k]
                    slength = genevertices[i][k+1]-genevertices[i][k]
                    sstart2 = start2 + (sstart1 - start1)
                    segment_matches_matrix_[i][j].append([sstart1,sstart2,slength,idty])
                    k += 1
            
    for i in range(len(genevertices)):
        for j in range(i,len(genevertices)):
            segment_matches_matrix_[i][j].sort()
            for k in range(len(segment_matches_matrix_[i][j])-1,0,-1):
                matchk = segment_matches_matrix_[i][j][k]
                matchkm1 = segment_matches_matrix_[i][j][k-1]
                if( matchk[:3] == matchkm1[:3]):
                    segment_matches_matrix_[i][j][k-1][3] = max(segment_matches_matrix_[i][j][k-1][3],segment_matches_matrix_[i][j][k][3])
                    del(segment_matches_matrix_[i][j][k])

    return segment_matches_matrix_

def compute_pairwise_scores(genevertices,genevertices2index,segment_matches_matrix,pscore):
    pairwise_score_matrix = []
    pairwise_score_matrix_ = []
    for i in range(len(genevertices)):
        pairwise_score_matrix.append([])
        pairwise_score_matrix_.append([])
        for j in range(len(genevertices)):
            pairwise_score_matrix[i].append([])
            pairwise_score_matrix_[i].append([])
    for i in range(len(genevertices)):
        for j in range(i,len(genevertices)):
            for k in range(len(genevertices2index[i])):
                pairwise_score_matrix[i][j].append([])
                pairwise_score_matrix_[i][j].append([])
                for l in range(len(genevertices2index[j])):
                    pairwise_score_matrix[i][j][k].append(0)
                    pairwise_score_matrix_[i][j][k].append(0)

    geneverticesimage = []
    for i in range(len(genevertices)):
        geneverticesimage.append([])
        for k in range(len(genevertices2index[i])):
            geneverticesimage[i].append([])

    for i in range(len(genevertices)):
        for j in range(i,len(genevertices)):
            for match in segment_matches_matrix[i][j]:
                start1,start2,length,idty= match
                sc =  0
                if(pscore == "unit"):
                    sc = 1
                else: #pscore == "idty"
                    sc = idty
                pairwise_score_matrix[i][j][genevertices2index[i][start1]][genevertices2index[j][start2]] = sc
                pairwise_score_matrix_[i][j][genevertices2index[i][start1]][genevertices2index[j][start2]] = sc
                geneverticesimage[i][genevertices2index[i][start1]].append((j,genevertices2index[j][start2]))
                geneverticesimage[j][genevertices2index[j][start2]].append((i,genevertices2index[i][start1]))
    
    for i in range(len(genevertices)):
        for j in range(i+1,len(genevertices)):
            for k in range(len(genevertices2index[i])):
                for l in range(len(genevertices2index[j])):
                    image_k =  geneverticesimage[i][k]
                    image_l = geneverticesimage[j][l]
                    inters = set(image_k).intersection(set(image_l))
                    for segment in inters:
                        geneidx = segment[0]
                        startidx = segment[1]
                        
                        scorei = -1
                        scorej = -1
                        if(geneidx != i and geneidx != j):
                            if(i < geneidx):
                                scorei = pairwise_score_matrix[i][geneidx][k][startidx]
                            else:
                                scorei = pairwise_score_matrix[geneidx][i][startidx][k]
                            if(j < geneidx):
                                scorej = pairwise_score_matrix[j][geneidx][l][startidx]
                            else:
                                scorej = pairwise_score_matrix[geneidx][j][startidx][l]
                            
                           
                            if(pscore == "unit"):
                                pairwise_score_matrix_[i][j][k][l] += 1
                            else:# pscore == "idty"
                                 pairwise_score_matrix_[i][j][k][l] += max(scorei,scorej)
                                
                                
    return pairwise_score_matrix_

 
def compute_geneid2index(targetdata):
    geneid2index = {}
    for i  in range(len(targetdata)):
        geneid,geneseq = targetdata[i]
        geneid2index[geneid]=i
    return geneid2index

def compute_genevertices2index(genevertices):
    genevertices2index = []
    for i in range(len(genevertices)):
        genevertices2index.append({})
        for j in range(len(genevertices[i])-1):
            genevertices2index[i][genevertices[i][j]]=j
    return genevertices2index        
    
def compute_msa_recursif(pairwise_score_matrix,node,geneid2index):
    columnlist = []
    
    #if node is a leaf
    if(node.name != None):
        geneid  = node.name
        i = geneid2index[geneid]
        for j in range(len(pairwise_score_matrix[i][i][0])):
            columnlist.append([[i,j]])
        
    #if node is an internal
    else:
        # compute msa at left child
        columnlistleft  = compute_msa_recursif(pairwise_score_matrix,node.descendants[0],geneid2index)
        # compute msa at right child
        columnlistright = compute_msa_recursif(pairwise_score_matrix,node.descendants[1],geneid2index)
        # merge left and right into a single msa
        columnlist = align(columnlistleft,columnlistright,pairwise_score_matrix)
    
    return columnlist

def align(columnlistleft,columnlistright,pairwise_score_matrix):
    columnlist = []
    aln_matrix = []
    for i in range(len(columnlistleft)+1):
        aln_matrix.append([])
        for j in range(len(columnlistright)+1):
            aln_matrix[i].append(0)

    for i in range(len(columnlistleft)):
        for j in range(len(columnlistright)):
            subs = aln_matrix[i-1][j-1] + score(columnlistleft[i],columnlistright[j],pairwise_score_matrix)
            dele = aln_matrix[i-1][j]
            inse =  aln_matrix[i][j-1]
            opt = max(max(subs,dele),inse)
            aln_matrix[i][j] = opt

    i = len(columnlistleft)-1
    j = len(columnlistright)-1
    while(i != -1 and j != -1):
        subs = aln_matrix[i-1][j-1] + score(columnlistleft[i],columnlistright[j],pairwise_score_matrix)
        dele = aln_matrix[i-1][j]
        inse =  aln_matrix[i][j-1]
        if(aln_matrix[i][j] == inse):
            columnlist = [columnlistright[j]] + columnlist
            j -=1
        elif(aln_matrix[i][j] == dele):
            columnlist = [columnlistleft[i]] + columnlist
            i -= 1
        else:
            columnlist = [columnlistleft[i]+columnlistright[j]] + columnlist
            i -= 1
            j -= 1
    while(i != -1):
        columnlist = [columnlistleft[i]] + columnlist
        i -= 1
    while(j != -1):
        columnlist = [columnlistright[j]] + columnlist
        j -= 1
                
    return columnlist
    
def score(columnlistleft,columnlistright,pairwise_score_matrix):
    score_ = 0
    for left in columnlistleft:
        i,k = left
        for right in columnlistright:
            j,l = right
            if(i < j):
                score_ += pairwise_score_matrix[i][j][k][l]
            else:
                score_ += pairwise_score_matrix[j][i][l][k]
    return score_

def column2mblocklist(columnlist,gene_list,genevertices):
    mblocklist = []
    for column in columnlist:
        mblock = {}#[]
        for segment in column:
            geneidx,posidx= segment
            #mblock = [[gene_list[geneidx],genevertices[geneidx][posidx],genevertices[geneidx][posidx+1]]] + mblock
            mblock[gene_list[geneidx]] = [genevertices[geneidx][posidx],genevertices[geneidx][posidx+1]]
        mblocklist.append(mblock)
    return mblocklist

def addsource2mblocklist(mblocklist,gene_list,gene2cds,cdsexon,cds2geneexon):
    mblocklist_ = []
    currentcdsexonidx = {}
    for geneid in gene_list:
        for cdsid in gene2cds[geneid]:
             currentcdsexonidx[cdsid] = 0
    for i in range(len(mblocklist)):
        mblock = {}
        nbgenes = len(mblocklist[i].keys())
        for geneid in mblocklist[i].keys():
            mblock[geneid] = mblocklist[i][geneid]
            genesegment = mblocklist[i][geneid]
            for cdsid in gene2cds[geneid]:
                if(currentcdsexonidx[cdsid] < len(cdsexon[cdsid])):
                    currentcdsexon = cdsexon[cdsid][currentcdsexonidx[cdsid]]
                    genecdsexon = cds2geneexon[cdsid][currentcdsexon[0]]
                    while(genesegment[0] >=  genecdsexon[1] and currentcdsexonidx[cdsid] < len(cdsexon[cdsid])):
                        currentcdsexonidx[cdsid] += 1
                        currentcdsexon = cdsexon[cdsid][currentcdsexonidx[cdsid]]
                        genecdsexon = cds2geneexon[cdsid][currentcdsexon[0]]
                    if(genecdsexon[0] <= genesegment[0] and genesegment[1] <= genecdsexon[1]):
                        mblock[cdsid] = [currentcdsexon[0]+genesegment[0]-genecdsexon[0],currentcdsexon[0]+genesegment[1]-genecdsexon[0]]
                        if(genesegment[1] == genecdsexon[1]):
                            currentcdsexonidx[cdsid] += 1
        if(len(mblock.keys()) > nbgenes):
            mblocklist_.append(mblock)
    return mblocklist_
    
def mergemblocks(mblocklist):
    i = len(mblocklist)-1
    while(i > 0):
        if(set(mblocklist[i].keys()) == set(mblocklist[i-1].keys())):
            if (all([mblocklist[i][geneid][0] == mblocklist[i-1][geneid][1] for geneid in mblocklist[i].keys()])):
                for geneid in mblocklist[i].keys():
                    mblocklist[i-1][geneid][1] = mblocklist[i][geneid][1]
                del(mblocklist[i])
        i -= 1
    return mblocklist
            

def compute_msa(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,cds2geneid,cds2geneexon,gene2cds,pscore,tree,outputprefix):

    geneid2index = compute_geneid2index(targetdata)
    
    segment_matches_matrix = comparisonresults2segmentmatches(targetdata, nbinitialsource, extendedsourcedata,comparisonresults,comparisonresults_idty,cdsexon,cds2geneexon,geneid2index)

    gene_list = [gene[0] for gene in targetdata]

    genevertices = compute_refine_boundaries(segment_matches_matrix,gene_list)
    
    genevertices2index = compute_genevertices2index(genevertices)

    segment_matches_matrix = refine_segment_matches_matrix(segment_matches_matrix,genevertices)

    pairwise_score_matrix = compute_pairwise_scores(genevertices,genevertices2index,segment_matches_matrix, pscore)

    t = newick.loads(tree)

    columnlist = compute_msa_recursif(pairwise_score_matrix,t[0],geneid2index)
    
    mblocklist = column2mblocklist(columnlist,gene_list,genevertices)

    mblocklist = addsource2mblocklist(mblocklist,gene_list,gene2cds,cdsexon,cds2geneexon)
    
    mblocklist = mergemblocks(mblocklist)

    return mblocklist


