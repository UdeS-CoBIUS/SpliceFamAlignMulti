#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compare.py`` **module description**:

This module is the module that launches the comparison for all pairs of source CDS and target gene.

.. moduleauthor:: Aida Ouangraoua

2019

"""

import networkx as nx
from collections import Counter
from Bio import pairwise2
import multiprocessing
from functools import partial
from contextlib import contextmanager
from multiprocessing import Pool
import time

MAX_EXTREMITY_DIFF = 30

# Create an initial graph
# Nodes : segements of gene and cds aligned in blocks
# Edges : alignments between segments
def  create_graph(targetdata, nbinitialsource, extendedsourcedata, comparisonresults,comparisonresults_idty):
    graphid2ensemblid = {}   # dictionnary graphid to ensemblid
    ensemblid2graphid = {}   # dictionnary graphid to ensemblid
    allcdsseq = {}
    mblocklist_graph = nx.Graph()
    nodeid = 0
    i = 0
    # for each block in each alignment
    for gene in targetdata:
        geneid,geneseq = gene
        graphid2ensemblid["g"+str(i)] = geneid
        ensemblid2graphid[geneid] = "g"+str(i)
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            graphid2ensemblid["c"+str(j)] = cdsid
            ensemblid2graphid[cdsid] = "c"+str(j)
            allcdsseq["c"+str(j)] = cdsseq
            #print(i,j,geneid,cdsid)
            null,blocklist,null,null,null = comparisonresults[i][j]
            for k in range(len(blocklist)):
                block = blocklist[k]
                gstart,gend = block[2:]
                cstart, cend = block[:2]
                # create a node for gene segment
                genenode = "g"+str(i)+"_"+str(gstart)+"_"+str(gend)
                # create a node for cds segment
                cdsnode = "c"+str(j)+"_"+str(cstart)+"_"+str(cend)
                # add an edge between node
                mblocklist_graph.add_node(genenode,
                                          id = "g"+str(i),
                                          start = gstart,
                                          end = gend)
                mblocklist_graph.add_node(cdsnode,
                                          id = "c"+str(j),
                                          start = cstart,
                                          end = cend)
                mblocklist_graph.add_edge(genenode,cdsnode, idty = comparisonresults_idty[i][j][k], connectivity = -1)
                
        i += 1
        
    return mblocklist_graph, graphid2ensemblid,ensemblid2graphid, allcdsseq

#  Connect that correspond to the same segments (equivalent nodes)
def connect_equivalent_nodes(extendedsourcedata,targetdata,mblocklist_graph,nbinitialsource,geneexon,cdsexon):
    # two segments are equivalent if they overlap and
    # their overlap is at least half of each of them and
    # ---! and the difference between starts is less than 30
    # ---! and the difference between ends is less than 30

    nodesperseq = {}

    for i in range(len(targetdata)):
        nodesperseq["g"+str(i)] = []
        
    for j in range(nbinitialsource):
        nodesperseq["c"+str(j)] = []
        
    for n in mblocklist_graph.nodes():
        word = n.split("_")
        pos = [int(x) for x in word[1:]]
        nodesperseq[word[0]].append([n]+pos)
        
    # Connect equivalent gene segments
    for id in nodesperseq.keys():
        nodes = nodesperseq[id]
        nbnodes = len(nodes)
        for k in range(nbnodes):
            u,startu,endu = nodes[k]
            for l in range(k+1,nbnodes):
                v,startv,endv = nodes[l]
                if((startv < endu) and (startu < endv)):
                    start = max(startu,startv)
                    end = min(endu,endv)
                    if (startu == startv and endu == endv):
#                    if (((end-start+1) > (1.0/2)*(endu-startu+1)) and
#                        ((end-start+1) > (1.0/2)*(endv-startv+1))):
                        mblocklist_graph.add_edge(u,v, idty = 1.0, connectivity = 100)
                        #print("ici",u,v)
                    #else:
                        #print("conflict",u,v)
    return mblocklist_graph

#  Compute edges connectivity : if the edge (a,b) is removed,
# connectivity = length of the shortest path between a and b
def weight_edges(mblocklist_graph):
    idty = nx.get_edge_attributes(mblocklist_graph,'idty')
    connec = nx.get_edge_attributes(mblocklist_graph,'connectivity')
    for (a,b) in list(mblocklist_graph.edges()):
        if(get_feature(a,b,connec) == -1):
            sim = get_feature(a,b,idty)
            weight = 0
            mblocklist_graph.remove_edge(a,b)
            if(nx.has_path(mblocklist_graph, a, b)):
                #path = sorted(nx.all_shortest_paths(mblocklist_graph, a, b))[0]
                path = nx.shortest_path(mblocklist_graph, a, b)
                #count number of block alignments in path (connectivity < 100)
                for k in range(len(path)-1):
                    if(get_feature(path[k],path[k+1],connec) < 100):
                        weight += 1
                # transform weight measure into a similary measure
                # the higher, the more connected 
                weight = 100 - weight
                assert(weight > 0)
            mblocklist_graph.add_edge(a,b, idty = sim, connectivity = weight)
    return mblocklist_graph
    
    
#  Remove erroneous edges (blocks) to disconnect conflicting nodes
def poolDisconnect(cc,mblocklist_graph, nbinitialsource):

    temps = time.time()
    #Discard cc with high number of nodes
    if(len(cc) >=  nbinitialsource * 10):
        #print(time.time()-temps, "not disconnect one cc", len(cc))
        return []

    # compute induced subgraph and weight edges
    cc = sorted(list(cc))
    cc_graph = nx.Graph(mblocklist_graph.subgraph(cc))
    cc_graph = weight_edges(cc_graph)

    idty = nx.get_edge_attributes(cc_graph,'idty')
    connec = nx.get_edge_attributes(cc_graph,'connectivity')
    # for any node u
    for i in range(len(cc)):
        u = cc[i]
        motu = u.split("_")
        idu = motu[0]
        startu,endu = [int(x) for x in motu[1:]]
        # for any node v
        for j in range(i+1,len(cc)):
            v = cc[j]
            motv = v.split("_")
            idv = motv[0]
            startv,endv = [int(x) for x in motv[1:]]
            # and u and v are conflicting nodes (same sequence but different locations)
            if ((idu == idv) and (endu <= startv or endv <= startu)):
                #print(idu,startu,endu,startv,endv)
                #while exists a path between u and v
                while(nx.has_path(cc_graph, u, v)):
                    #find a shortest path
                    path = nx.shortest_path(cc_graph, u, v)
                    #print("...",path)
                    ## remove the smallest similarity edge among
                    ## the smallest connectivity edges on the path
                    min_connec = 100
                    min_sim = 1.0
                    sim = 1.0
                    a = path[0]
                    b = path[1]
                    for k in range(len(path)-1):
                        cur_connec =  get_feature(path[k],path[k+1],connec)
                        cur_sim =  get_feature(path[k],path[k+1],idty)
                        if(cur_connec < min_connec or (cur_connec == min_connec and cur_sim < min_sim)):
                            min_connec  = cur_connec
                            min_sim = cur_sim
                            a = path[k]
                            b = path[k+1]
                        if(cur_sim < sim):
                            sim = cur_sim
                    #print("...",a,b, min_connec,min_sim)                    
                    cc_graph.remove_edge(a,b)
    connectedcomponents = list(nx.connected_components(cc_graph))
    #print(time.time()-temps, "disconnect one cc", len(cc))
    return connectedcomponents

# return the value of feature "feature_name" for an edge
def get_feature(a,b,feature_name):
    if((a,b) in feature_name.keys()):
        feat =  feature_name[(a,b)]
    else:
        feat =  feature_name[(b,a)]
    return feat

# Compute initial mblocklist : 1 connected component (cc) = 1 mblock
# nodes from the same sequence in the same cc are merges into one entry in mblock
def pool_initialize_mblocklist(cc, graphid2ensemblid,cdsexon,geneexon):
    #temps = time.time()
    mblock = {}
    for u in cc:
        motu = u.split("_")
        idu = motu[0]
        startu,endu = [int(x) for x in motu[1:]]
        # if sequence already has an entry in current mblock
        if(idu in mblock.keys()):
            cur_startu,cur_endu = mblock[idu]
            # if current location is an exon location
            if((idu[0]=='c' and [cur_startu,cur_endu] in cdsexon[graphid2ensemblid[idu]])
               or (idu[0]=='g' and [cur_startu,cur_endu] in geneexon[graphid2ensemblid[idu]])):
                # if new location is also an exon location
                if((idu[0]=='c' and [startu,endu] in cdsexon[graphid2ensemblid[idu]])
                   or (idu[0]=='g' and [startu,endu] in geneexon[graphid2ensemblid[idu]])):
                    # keep the largest length location 
                    if(endu-startu > cur_endu-cur_startu):
                         mblock[idu] = [startu,endu]
            # else current location is not an exon location
            else:
                # if new location is an exon location
                if((idu[0]=='c' and [startu,endu] in cdsexon[graphid2ensemblid[idu]])
                   or (idu[0]=='g' and [startu,endu] in geneexon[graphid2ensemblid[idu]])):
                    # keep exon location
                    mblock[idu] = [startu,endu]
                # new location is not an exon location neither
                else:
                    # compute maximum length
                    mblock[idu][0]=min(startu,cur_startu)
                    mblock[idu][1]=max(endu,cur_endu)
        else:
            mblock[idu] = [startu,endu]
    #print(time.time() - temps)
    return mblock

# Sort mblock by decreasing number of entries
def sort_mblocklist(mblocklist):
    for i in range(len(mblocklist)):
        for j in range(1,len(mblocklist) - i):
            if(len(mblocklist[j-1].keys()) < len(mblocklist[j].keys())):
                tmp = mblocklist[j-1]
                mblocklist[j-1] = mblocklist[j]
                mblocklist[j] = tmp
            #for determinism
            elif(len(mblocklist[j-1].keys()) == len(mblocklist[j].keys())):
                id1 = "".join(sorted(list(mblocklist[j-1].keys())))
                id2 = "".join(sorted(list(mblocklist[j].keys())))
                if(id1 < id2):
                    tmp = mblocklist[j-1]
                    mblocklist[j-1] = mblocklist[j]
                    mblocklist[j] = tmp
                elif(id1 == id2):
                    id = sorted(list(mblocklist[j-1].keys()))[0]
                    if(mblocklist[j-1][id][0] <  mblocklist[j][id][0]):
                        tmp = mblocklist[j-1]
                        mblocklist[j-1] = mblocklist[j]
                        mblocklist[j] = tmp
                        
    return mblocklist
    

# Create mblocklist by progressively adding mblocks (decreasing number of entries)
def create_mblocklist(mblocklist_init,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid,cdsexon,geneexon):
    nbremoved = 0
    mblocklist = []
    # Add mblocks supported by at least 3 entries
    for i in range(len(mblocklist_init)):
        if(len(mblocklist_init[i].keys()) > 2):
            mbocklist, status, initial_nbid, final_nbid  = add_mblock(mblocklist,mblocklist_init[i],graphid2ensemblid,cdsexon,geneexon)

            mblocklist = order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)
            if(status == "remove"):
                #print(mblocklist_init[i])
                nbremoved += 1
            #else:
                #print(initial_nbid, final_nbid)
    return mblocklist, nbremoved

# Add mblock in mblocklist
def add_mblock(mblocklist,mblock,graphid2ensemblid,cdsexon,geneexon):
    initial_nbid = len(mblock.keys())
    status = ""
    if(len(mblocklist) == 0):
        mblocklist.append(mblock)
    else:
        position = []
        # compute locations of mblock according to already added mblocks 
        for i in range(len(mblocklist)):
            before,after,included,contain,overlap_before, overlap_after,overlap_ratio = compare_location(mblock,mblocklist[i])
            position.append([before,after,included,contain,overlap_before, overlap_after,overlap_ratio])

        mblock, pos, status = find_location(mblock,mblocklist,position,graphid2ensemblid,cdsexon,geneexon)

        # if location is before some location pos
        if(status == "before"):
            # remove all entries in mblock that contradicts this location
            for i in range(pos,len(mblocklist)):
                for k in range(1,6):
                    for id in position[i][k]:
                        if(id in mblock.keys() and  mblock[id][1] >  mblocklist[i][id][0]):
                            mblock.pop(id)
            #insert mblock at location
            if(len(mblock.keys()) > 1):
                mblocklist.insert(pos,mblock)
        # if position is at the end of mblocklist
        elif(status == "after"):
            #insert mblock at the end
            if(len(mblock.keys()) > 1):
                mblocklist.append(mblock)
        # if position is included or contains another mblock at location pos
        elif(status == "included" or status == "contain"):
            # remove all entries in mblock that contradicts this location
            for i in range(pos+1,len(mblocklist)):
                for k in range(1,6):
                    for id in position[i][k]:
                        if(id in mblock.keys() and  mblock[id][1] >  mblocklist[i][id][0]):
                            mblock.pop(id)
            #merge the two mblocks
            if(len(mblock.keys()) > 1):
                for id in mblock.keys():
                    if(id in mblocklist[pos].keys()):
                        mblocklist[pos][id][0] = min(mblocklist[pos][id][0], mblock[id][0])
                        mblocklist[pos][id][1] = max(mblocklist[pos][id][1], mblock[id][1])
                    else:
                        mblocklist[pos][id] = mblock[id]

    final_nbid = len(mblock.keys())
    return mblocklist, status, initial_nbid, final_nbid

# Find location of a new  mblock in mblocklist given the positions
#relative to each mblock in the mblocklist
def find_location(mblock,mblocklist,position,graphid2ensemblid,cdsexon,geneexon):
    status = "after"
    i = -1
    # go through mblocks while position has been found
    # and for each mblock do:
    while (status == "after" and i < (len(mblocklist)-1)):
        i += 1
        #comute number of common id per type of location
        nb_per_status = [len(position[i][k]) for k in range(6)]
        
        # compute total number of common id
        nb_total = sum(nb_per_status)

        # if no common id, continue to next mblock
        if(nb_total == 0):
            continue
        
        #compute ratio of id per type of location
        ratio_per_status= [1.0*nb_per_status[k]/nb_total for k in range(6)]

        ratio_before, ratio_after, ratio_included, ratio_contain, ratio_overlap_before, ratio_overlap_after = ratio_per_status

        #compute number of types of location with at least one id
        nb_nonnul = int(sum([nb_per_status[k]*1/(nb_per_status[k]-0.001) for k in range(6)]))

        #comute average rate of overlap before
        overlap_b_rate = 0.0
        nb_id = 0
        for id in position[i][4]:
            if(id in mblock.keys()):
                nb_id += 1 
                overlap_b_rate += position[i][6][id][0]
                #1.0*(mblock[id][1]-mblocklist[i][id][0])/(mblock[id][1]-mblock[id][0])
        if(nb_id > 0):
            overlap_b_rate = overlap_b_rate/nb_id
            
        #compute min rate of overlap after
        overlap_a_rate = 0.0
        nb_id = 0
        for id in position[i][5]:
            if(id in mblock.keys()):
                nb_id += 1
                overlap_a_rate += position[i][6][id][0]
                #1.0*(mblocklist[i][id][1]-mblock[id][0])/(mblock[id][1]-mblock[id][0])
        if(nb_id > 0):
            overlap_a_rate = overlap_a_rate/nb_id

        #compute min the rate of real exon replaced by a larger block for the location type "contain"
        erase_exon_rate = 0.0
        nb_id = 0
        for id in position[i][3]:
            if(id in mblock.keys()):
                nb_id += 1
                if (((id[0]=='c') and (not(mblock[id] in cdsexon[graphid2ensemblid[id]])) and (mblocklist[i][id] in cdsexon[graphid2ensemblid[id]]))
                    or ((id[0]=='g') and (not(mblock[id] in geneexon[graphid2ensemblid[id]])) and (mblocklist[i][id] in geneexon[graphid2ensemblid[id]]))):
                    erase_exon_rate += 1
        if(nb_id > 0):
            erase_exon_rate = erase_exon_rate/nb_id

        # if more than two possible locations, do not return any location
        if(nb_nonnul > 2):
            status = "remove"
        # if location is before or overlap_before with overlap_rate <= 0.25 : merge both cases in before
        elif(ratio_before == max(ratio_per_status) or (ratio_overlap_before == max(ratio_per_status) and overlap_b_rate <= 0.25)):
            # remove all entries in mblock that contradicts this location
            for k in list(range(1,4))+list(range(5,6)):
                for id in position[i][k]:
                    if(id in mblock.keys()):
                        mblock.pop(id)
            # correct all entries with overlapping locations to remove overlap
            for id in position[i][4]:
                if(id in mblock.keys()):
                    mblock[id] = [mblock[id][0],mblocklist[i][id][0]]
            ratio_before += ratio_overlap_before
            ratio_overlap_before = 0.0            
            status = "before"
            
        # if location is after or overlap_after with overlap_rate <= 0.25 : symmetric to previous case
        elif(ratio_after == max(ratio_per_status)  or (ratio_overlap_after == max(ratio_per_status) and overlap_a_rate < 0.25)):
            for k in list(range(1)) + list(range(2,5)):
                for id in position[i][k]:
                    if(id in mblock.keys()):
                        mblock.pop(id)
            for id in position[i][5]:
                if(id in mblock.keys()):
                    mblock[id] = [mblocklist[i][id][1],mblock[id][1]]
            ratio_after += ratio_overlap_after
            ratio_overlap_after = 0.0
            status = "after"
            
        # if location is include 
        elif(ratio_included == max(ratio_per_status)):
            # remove all entries in mblock that contradicts this location
            for k in list(range(2)) + list(range(3,6)):
                for id in position[i][k]:
                    if(id in mblock.keys()):
                        mblock.pop(id)
            status = "included"
            
        # if location is contain and no exon erased
        elif(ratio_contain == max(ratio_per_status) and erase_exon_rate == 0.0):
            # remove all entries in mblock that contradicts this location
            for k in list(range(3)) + list(range(4,6)):
                for id in position[i][k]:
                    if(id in mblock.keys()):
                        mblock.pop(id)
            status = "contain"
        #else do not return any location

        else:                        
            status = "remove"
        #print(ratio_before, ratio_after, ratio_included, ratio_contain, ratio_overlap_before, ratio_overlap_after)
            
    return mblock, i, status
        
# Compute the location of a mblock according to a mblock_i
def compare_location(mblock,mblock_i):
    before,after,included,contain,overlap_before,overlap_after,overlap_ratio = [],[],[],[],[],[],{}
    common_keys = sorted(list(set(mblock.keys()).intersection(mblock_i.keys())))
    for key in common_keys:
        mblock_length = mblock[key][1]-mblock[key][0]
        mblock_i_length = mblock_i[key][1]-mblock_i[key][0]

        if(mblock[key][1] <= mblock_i[key][0]):
            before.append(key)
        elif(mblock_i[key][1] <= mblock[key][0]):
            after.append(key)
        elif(mblock_i[key][0] <= mblock[key][0] and mblock[key][1] <= mblock_i[key][1]):
            included.append(key)
        elif(mblock[key][0] <= mblock_i[key][0] and mblock_i[key][1] <= mblock[key][1]):
            contain.append(key)
        elif(mblock[key][0] <= mblock_i[key][0]):
            overlap_before.append(key)
            overlap_length = mblock[key][1]-mblock_i[key][0]
            overlap_ratio[key]=[1.0*overlap_length/mblock_length,1.0*overlap_length/mblock_i_length]
        else:
            overlap_after.append(key)
            overlap_length = mblock_i[key][1]-mblock[key][0]
            overlap_ratio[key]=[1.0*overlap_length/mblock_length,1.0*overlap_length/mblock_i_length]
    return before,after,included,contain,overlap_before, overlap_after,overlap_ratio

# Remove mblock with only only gene entries
def remove_geneandsinglemblocks(mblocklist):
    to_delete = []
    for i in range(len(mblocklist)):
        #if(len(mblocklist[i].keys()) <= 1):
        #    to_delete.append(i)
        if(all([id[0]=='g' for id in mblocklist[i].keys()])):
            to_delete.append(i)

    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])
    return mblocklist

# Add cds region that are not present in alignment as new or extended mblocks 
def  complete_cds(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,cds2geneid,cds2geneexon):
    for cdsid in cdsexon.keys():
        graphid = ensemblid2graphid[cdsid]
        i = 0
        for exon in cdsexon[cdsid]:
            while(i < len(mblocklist) and ((not(graphid in mblocklist[i].keys())) or (mblocklist[i][graphid][1] <= exon[0]))):
                i += 1
            # if cds exon not in alignment until end: add new mblock
            if(i == len(mblocklist)):
                #print("add_at_end",cdsid,exon)
                mblock = {}
                mblock[graphid] = exon
                geneid = cds2geneid[cdsid]
                graphgeneid = ensemblid2graphid[geneid]
                cdstart = exon[0]
                #mblock[graphgeneid] = cds2geneexon[cdsid][cdstart]
                mblocklist.append(mblock)
            # if cds exon in alignment : continue
            elif(exon == mblocklist[i][graphid]):
                continue
            # if cds exon not in alignment until position i that is after: add new mblock
            elif(exon[1] <= mblocklist[i][graphid][0]):
                #print("add_before_i",cdsid,exon)
                mblock = {}
                mblock[graphid] = exon
                geneid = cds2geneid[cdsid]
                graphgeneid = ensemblid2graphid[geneid]
                cdstart = exon[0]
                #mblock[graphgeneid] = cds2geneexon[cdsid][cdstart]
                mblocklist.insert(i,mblock)
            # if overlap region in alignment : extend to cover exon
            else:
                #print("replace_i",cdsid,exon,mblocklist[i][graphid])
                geneid = cds2geneid[cdsid]
                graphgeneid = ensemblid2graphid[geneid]
                cdstart = exon[0]
                exongene = cds2geneexon[cdsid][cdstart]
                if (not(graphgeneid in list(mblocklist[i].keys()))):
                    #mblocklist[i][graphgeneid] = exongene
                    cdstart = exon[0]
                if(exon[0] < mblocklist[i][graphid][0]):
                     mblocklist[i][graphid][0] = exon[0]
                     #mblocklist[i][graphgeneid][0] = exongene[0]
                if(mblocklist[i][graphid][1] < exon[1]):
                     mblocklist[i][graphid][1] = exon[1]
                     #mblocklist[i][graphgeneid][1] = exongene[1]
                     
    return mblocklist

# Add or correct gene location corresponding to cds location
def add_genelocation(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,geneexon,cds2geneid,cds2geneexon):
    for i in range(len(mblocklist)):
        gene_location = {}
        cdsperlocation = {}

        for graphid in mblocklist[i].keys():
            if(graphid[0]=='c'):
                cdsid = graphid2ensemblid[graphid]
                geneid = cds2geneid[cdsid]
                graphgeneid = ensemblid2graphid[geneid]
                if(graphgeneid not in gene_location.keys()):
                    gene_location[graphgeneid] = []
                    cdsperlocation[graphgeneid] = {}
                gene_location[graphgeneid].append(cds2genelocation(cdsid,mblocklist[i][graphid],cdsexon,cds2geneexon))
                gstart,gend = gene_location[graphgeneid][-1]
                locationid = str(gstart)+"_"+str(gend)
                if(locationid not in cdsperlocation[graphgeneid].keys()):
                    cdsperlocation[graphgeneid][locationid] = []
                cdsperlocation[graphgeneid][locationid].append(graphid)
                    
        for graphgeneid in gene_location.keys():
            maxcds = 0
            genelocation = [0,0]
            for pos in cdsperlocation[graphgeneid].keys():
                if (len(cdsperlocation[graphgeneid][pos]) > maxcds):
                    maxcds  = len(cdsperlocation[graphgeneid][pos])
                    genelocation = pos
            if(genelocation.split("_")[0] != "-1" and genelocation.split("_")[1] != "-1"):
                [gstart,gend] = [int(x) for x in genelocation.split("_")]
                mblocklist[i][graphgeneid] = [gstart,gend]
                for pos in cdsperlocation[graphgeneid].keys():
                    if(pos != genelocation):
                        [pstart,pend] = [int(x) for x in pos.split("_")]
                        for graphid in cdsperlocation[graphgeneid][pos]:
                            # locations are overlapping, and current location is not an exon
                            if(pstart < gend and gstart < pend):
                                if([pstart,pend] not in geneexon[graphid2ensemblid[graphgeneid]]):
                                    cstart = mblocklist[i][graphid][0] + (gstart - pstart)
                                    cend = mblocklist[i][graphid][1] + (gend - pend)
                                    mblocklist[i][graphid] = [cstart,pend]
                            else:
                                mblocklist[i].pop(graphid)
    return mblocklist


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
    
    
# Merge overlapping gene blocks
def merge_overlapping(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid,cdsexon,geneexon):
    to_delete = []
    # For any pair of mblocks
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            if(i in to_delete or j in to_delete):
                continue
            overlap = []
            disjoint = []
            # compute sets of overlap/disjoint ids
            for id in mblocklist[i].keys():
                if(id in mblocklist[j].keys()):
                    if (not((mblocklist[i][id][1] <= mblocklist[j][id][0]) or (mblocklist[j][id][1] <= mblocklist[i][id][0]))):
                        overlap.append(id)
                    else:
                        disjoint.append(id)
            # if there are conflicts
            if(len(overlap) > 0 and len(disjoint) > 0):
                #remove gene without cds and recompute sets of overlap/disjoint
                overlap = []
                disjoint = []
                for id in overlap + disjoint:
                    if(id[0] == "g"):
                        keepj = False
                        keepi = False
                        for cid in [ensemblid2graphid[x] for x in gene2cds[graphid2ensemblid[id]]]:
                            if((cid in mblocklist[j].keys())):
                                keepj = True
                                    
                            if((cid in mblocklist[i].keys())):
                                keepi = True

                        if(not keepj):
                            mblocklist[j].pop(id)
                        if(not keepi):
                            mblocklist[i].pop(id)
                            
                        if ((id in mblocklist[i].keys()) and (id in mblocklist[j].keys()) and not((mblocklist[i][id][1] <= mblocklist[j][id][0]) or (mblocklist[j][id][1] <= mblocklist[i][id][0]))):
                            overlap.append(id)
                        else:
                            disjoint.append(id)

            # if there are still conflicts
            if(len(overlap) > 0 and len(disjoint) > 0):
                switch = []
                # merge close disjoint segments
                for id in disjoint:
                    if(id[0] == 'g' and (abs(mblocklist[i][id][1]- mblocklist[j][id][0]) <= 3 or abs(mblocklist[i][id][0]- mblocklist[j][id][1]) <= 3)):
                        switch.append(id)
                    elif(id[0] == 'c'  and (mblocklist[i][id][1] ==  mblocklist[j][id][0] or mblocklist[i][id][0] == mblocklist[j][id][1])):
                        switch.append(id)
                for id in switch:
                    overlap.append(id)
                    disjoint.remove(id)

            if(len(overlap) > 0):
                # if there is no conflict, merge mblocks
                if(len(disjoint) == 0):

                    for id in mblocklist[j].keys():
                        if(id in mblocklist[i].keys()):
                            mblocklist[i][id]=[min(mblocklist[i][id][0],mblocklist[j][id][0]),max(mblocklist[i][id][1],mblocklist[j][id][1])]
                        else:
                             mblocklist[i][id]=mblocklist[j][id]
                    to_delete.append(j)
                # otherwise
                else:
                    # delete the mblock with the smallest number of exons
                    nb_exon_i = 0
                    nb_exon_j = 0
                    #print("delete", overlap, disjoint, "\n", mblocklist[i], mblocklist[j])
                    for id in overlap+disjoint:
                        if id[0]=='g':
                            geneid = graphid2ensemblid[id]
                            if(mblocklist[i][id] in geneexon[geneid]):
                                nb_exon_i += 1
                            if(mblocklist[j][id] in geneexon[geneid]):
                                nb_exon_j += 1
                        else:
                            cdsid = graphid2ensemblid[id]
                            if(mblocklist[i][id] in cdsexon[cdsid]):
                                nb_exon_i += 1
                            if(mblocklist[j][id] in cdsexon[cdsid]):
                                nb_exon_j += 1

                    if(nb_exon_i > nb_exon_j):
                        to_delete.append(j)
                    elif(nb_exon_j > nb_exon_i):
                        #to_delete.append(i)
                        tmp = mblocklist[j]
                        mblocklist[j] = mblocklist[i]
                        mblocklist[i] = tmp
                        to_delete.append(j)
                    else:
                        id = overlap[0]
                        length_i = mblocklist[i][id][1] - mblocklist[i][id][0]
                        length_j = mblocklist[j][id][1] - mblocklist[j][id][0]
                        if(length_i > length_j):
                            to_delete.append(j)
                        else:
                            #to_delete.append(i)
                            tmp = mblocklist[j]
                            mblocklist[j] = mblocklist[i]
                            mblocklist[i] = tmp
                            to_delete.append(j)
                    #print(mblocklist[j])
                                        
    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])

                        
    return mblocklist

# Replace graphid by ensemblid
def  changesequenceid(mblocklist,graphid2ensemblid):
    for i in range(len(mblocklist)):
        mblock = {}
        for graphid in mblocklist[i].keys():
            mblock[graphid2ensemblid[graphid]] =  mblocklist[i][graphid]
        del(mblocklist[i])
        mblocklist.insert(i,mblock)

    return mblocklist

def order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid):
    """
    This function orders blocks in the blocklist by increasing
    order of query start location.
    """
    
    for i in range(len(mblocklist)):
        min = i
        for j in range(i+1,len(mblocklist)):
            inf = []
            sup = []
            eq = []
            common_keys = list(set(mblocklist[j].keys()).intersection(mblocklist[min].keys()))
            
            for id in common_keys:
                if(mblocklist[j][id][0] <  mblocklist[min][id][0]):
                    inf.append(id)
                elif(mblocklist[min][id][0] <  mblocklist[j][id][0]):
                    sup.append(id)
                else:
                    eq.append(id)
            if(len(inf) > len(sup + eq)):
                # if(len(sup + eq) >  0):
                    # print("inf",len(sup), len(inf), len(eq),"conflicts in order")
                    # print(inf[0],mblocklist[j][inf[0]][0], mblocklist[min][inf[0]][0])
                    # if(len(sup) >  0):
                    #     print(sup[0],mblocklist[j][sup[0]][0], mblocklist[min][sup[0]][0])
                    # if(len(eq) >  0):
                    #     print(seq[0],mblocklist[j][eq[0]][0], mblocklist[min][eq[0]][0])
                for id in sup + eq:
                    if(id[0] == "g"):
                        keepj = False
                        keepmin = False
                        for cid in [ensemblid2graphid[x] for x in gene2cds[graphid2ensemblid[id]]]:
                            if((cid in mblocklist[j].keys())):
                                keepj = True
                                
                            if((cid in mblocklist[min].keys())):
                                keepmin = True

                        if(not keepj):
                            # print("pop", id, mblocklist[j][id])
                            mblocklist[j].pop(id)
                        if(not keepmin):
                            # print("pop", id, mblocklist[min][id])
                            mblocklist[min].pop(id)
                    else:
                        # print("pop", id, mblocklist[j][id])
                        # print("pop", id, mblocklist[min][id])
                        mblocklist[j].pop(id)
                        mblocklist[min].pop(id)
                min = j
            elif(len(sup) > len(inf + eq)):
                # if(len(inf + eq) >  0):
                    # print("sup", len(sup), len(inf), len(eq), "conflicts in order")
                    # print(sup[0],mblocklist[j][sup[0]][0], mblocklist[min][sup[0]][0])
                    # if(len(inf) >  0):
                    #     print(inf[0],mblocklist[j][inf[0]][0], mblocklist[min][inf[0]][0])
                    # if(len(eq) >  0):
                    #     print(eq[0],mblocklist[j][eq[0]][0], mblocklist[min][eq[0]][0])
                for id in inf + eq:
                    if(id[0] == "g"):
                        keepj = False
                        keepmin = False
                        for cid in [ensemblid2graphid[x] for x in gene2cds[graphid2ensemblid[id]]]:
                            if((cid in mblocklist[j].keys())):
                                keepj = True
                                
                            if((cid in mblocklist[min].keys())):
                                keepmin = True

                        if(not keepj):
                            # print("pop", id, mblocklist[j][id])
                            mblocklist[j].pop(id)
                        if(not keepmin):
                            # print("pop", id, mblocklist[min][id])
                            mblocklist[min].pop(id)
                    else:
                        # print("pop", id, mblocklist[j][id])
                        # print("pop", id, mblocklist[min][id])
                        mblocklist[j].pop(id)
                        mblocklist[min].pop(id)
                        
            elif(len(eq) > len(inf + sup)):
                # if(len(inf + sup) >  0):
                    # print("eq", len(sup), len(inf), len(eq),"conflicts in order")
                    # print(eq[0],mblocklist[j][eq[0]][0], mblocklist[min][eq[0]][0])
                    # if(len(sup) >  0):
                    #     print(sup[0],mblocklist[j][sup[0]][0], mblocklist[min][sup[0]][0])
                    # if(len(inf) >  0):
                    #     print(inf[0],mblocklist[j][inf[0]][0], mblocklist[min][inf[0]][0])
                for id in inf + sup:
                    if(id[0] == "g"):
                        keepj = False
                        keepmin = False
                        for cid in [ensemblid2graphid[x] for x in gene2cds[graphid2ensemblid[id]]]:
                            if((cid in mblocklist[j].keys())):
                                keepj = True
                                
                            if((cid in mblocklist[min].keys())):
                                keepmin = True

                        if(not keepj):
                            # print("pop", id, mblocklist[j][id])
                            mblocklist[j].pop(id)
                        if(not keepmin):
                            # print("pop",id,  mblocklist[min][id])
                            mblocklist[min].pop(id)
                    else:
                        # print("pop", id, mblocklist[j][id])
                        # print("pop", id, mblocklist[min][id])
                        mblocklist[j].pop(id)
                        mblocklist[min].pop(id)
                
            else:
                # if(len(sup+inf+eq)!= 0):
                    # print(len(sup), len(inf), len(eq), "no decision in order")
                for id in inf+sup+eq:
                    if(id[0] == "g"):
                        keepj = False
                        keepmin = False
                        for cid in [ensemblid2graphid[x] for x in gene2cds[graphid2ensemblid[id]]]:
                            if((cid in mblocklist[j].keys())):
                                keepj = True
                                
                            if((cid in mblocklist[min].keys())):
                                keepmin = True

                        if(not keepj):
                            # print("pop", id, mblocklist[j][id])
                            mblocklist[j].pop(id)
                        if(not keepmin):
                            # print("pop", id, mblocklist[min][id])
                            mblocklist[min].pop(id)
                    else:
                        # print("pop", id, mblocklist[j][id])
                        # print("pop", id, mblocklist[min][id])
                        mblocklist[j].pop(id)
                        mblocklist[min].pop(id)
        if(min != i):
            tmp = mblocklist[min]
            mblocklist[min] = mblocklist[i]
            mblocklist[i] = tmp
    return mblocklist

def complete_and_merge(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,geneexon,cds2geneid,cds2geneexon,gene2cds):

    mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)

    
    mblocklist = complete_cds(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,cds2geneid,cds2geneexon)


    mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)


    mblocklist = add_genelocation(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,geneexon,cds2geneid,cds2geneexon)

    mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)


    mblocklist = merge_overlapping(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid,cdsexon,geneexon)


    mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)


    mblocklist = remove_geneandsinglemblocks(mblocklist)

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
    for k in range(10):
        for i in range(len(mblocklist)):
            for j in (order[i][before]-set([-1,len(mblocklist)])):
                order[i][before] = (order[i][before])|(order[j][before])
            for j in (order[i][after]-set([-1,len(mblocklist)])):
                order[i][after] = (order[i][after])|(order[j][after])
    return order

def mblocklength(mblock,allcdsseq):
    lengths = []
    for id in mblock.keys():
        if(id[0]=='c'):
            lengths.append(mblock[id][1]-mblock[id][0])
    length,count = Counter(lengths).most_common(1)[0]
    for id in mblock.keys():
        if(id[0]=='c' and (mblock[id][1]-mblock[id][0]) == length):
            seq = allcdsseq[id][mblock[id][0]:mblock[id][1]]
            break
    return length,seq

def cdslist(mblock):
    cds = []
    for id in mblock.keys():
        if(id[0]=='c'):
            cds.append(id)
    return cds

def merge_compatible_unordered(mblocklist,allcdsseq):
    mblocklist = remove_geneandsinglemblocks(mblocklist)

    before = 0
    after = 1
    order = partial_order(mblocklist)
    to_delete = []
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            if(i in to_delete or j in to_delete):
                continue
            if(not(j in (order[i][before])|(order[i][after]))):
                pos_before = max(max(order[i][before]),max(order[j][before]))
                pos_after = min(min(order[i][after]),min(order[j][after]))
                
                if (pos_before < pos_after):
                    lengthi,seqi = mblocklength(mblocklist[i],allcdsseq)
                    lengthj,seqj = mblocklength(mblocklist[j],allcdsseq)
                    if(abs(lengthi - lengthj) < 7 and lengthi%3 == lengthj%3 and compute_sequence_identity(seqi,seqj) > .3):
                        for id in mblocklist[j].keys():
                            mblocklist[i][id]=mblocklist[j][id]
                        order[i][before] = order[i][before]|order[j][before]
                        order[i][after] = order[i][after]|order[j][after]
                        to_delete.append(j)

    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])
    
    return mblocklist

def merge_compatible_ordered(mblocklist,allcdsseq):
    mblocklist = remove_geneandsinglemblocks(mblocklist)
        
    before = 0
    after = 1
    order = partial_order(mblocklist)
    to_delete = []
    for i in range(len(mblocklist)-1):
        if(i in to_delete):
            continue
        lengthi,seqi = mblocklength(mblocklist[i],allcdsseq)

        compatiblefound = False
        j = i+1
        
        while((not compatiblefound) and j < len(mblocklist)):            
            if(j in to_delete):
                j += 1
                continue
            lengthj,seqj = mblocklength(mblocklist[j],allcdsseq)

            if((min(order[i][after]) == j or max(order[j][before]) == i) and len(set(cdslist(mblocklist[i])) & set(cdslist(mblocklist[j]))) == 0 and abs(lengthi - lengthj) < 7 and lengthi%3 == lengthj%3 and compute_sequence_identity(seqi,seqj) > .3):
                for id in mblocklist[j].keys():
                    mblocklist[i][id]=mblocklist[j][id]
                to_delete.append(j)
                compatiblefound = True
            j += 1
            
    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])

    return mblocklist

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


def merge_compatible_extremity(mblocklist):
    mblocklist = remove_geneandsinglemblocks(mblocklist)

    before = 0
    after = 1
    order = partial_order(mblocklist)
    to_delete = []
    for i in range(len(mblocklist)):
        for j in range(i+1,len(mblocklist)):
            if(i in to_delete or j in to_delete):
                continue
            if(not(j in (order[i][before])|(order[i][after]))):
                if((order[i][before]== order[j][before] == set([-1]))  or
                   (order[i][after]== order[j][after] == set([len(mblocklist)]))):
                    #lengthi = mblocklength(mblocklist[i])
                    #lengthj = mblocklength(mblocklist[j])
                    #if(lengthi%3 == lengthj%3):
                    for id in mblocklist[j].keys():
                        mblocklist[i][id]=mblocklist[j][id]
                    to_delete.append(j)

    to_delete.sort(reverse = True)
    for i in to_delete:
        del(mblocklist[i])

    return mblocklist

                
def compute_msa(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,cds2geneid,cds2geneexon,gene2cds, compareExon):

    temps=time.time()

    mblocklist_graph, graphid2ensemblid,ensemblid2graphid,allcdsseq = create_graph(targetdata, nbinitialsource, extendedsourcedata, comparisonresults,comparisonresults_idty)

    #connected_components = list(nx.connected_components(mblocklist_graph))
    #print(time.time()-temps, "\n" ,len(connected_components), "initial connected components\n")

    mblocklist_graph = connect_equivalent_nodes(extendedsourcedata,targetdata,mblocklist_graph,nbinitialsource,geneexon,cdsexon)

    connected_components = list(nx.connected_components(mblocklist_graph))
    #print(time.time()-temps, "\n" ,len(connected_components), "connected components after connection of equivalent nodes\n")

    new_connected_components = []
    for cc in connected_components:
        nodesperseq = {}
        for n in list(cc):
            word = n.split("_")
            nodesperseq[word[0]] = 1
        if(len(cc) > 2 and len(cc) > len(nodesperseq.keys())):
            new_connected_components += poolDisconnect(cc,mblocklist_graph,nbinitialsource)
        else:
            new_connected_components.append(cc)
    connected_components = new_connected_components

    print(time.time()-temps, "\n" ,len(connected_components), "connected components after disconnection of conflicts\n")

    mblocklist_init = []
    for cc in connected_components:
        mblocklist_init.append(pool_initialize_mblocklist(cc, graphid2ensemblid,cdsexon,geneexon))
    mblocklist_init = sort_mblocklist(mblocklist_init) 


    #print(time.time()-temps, "\n" ,len(mblocklist_init), "initial mblocks\n")

    mblocklist_init = remove_geneandsinglemblocks(mblocklist_init)

    #print(time.time()-temps, "\n" ,len(mblocklist_init), "initial mblocks after removing gene mblocks\n")

    mblocklist,nbremoved = create_mblocklist(mblocklist_init,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid,cdsexon,geneexon)
    #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after create mblocks,", nbremoved, "initial mblocks removed\n")
    
    mblocklist = remove_geneandsinglemblocks(mblocklist)

    #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after removing gene mblocks\n")

    mblocklist = complete_and_merge(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,geneexon,cds2geneid,cds2geneexon,gene2cds)

    #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after completing cds, merging overlaps, and removing gene mblocks\n")

    mblocklist = complete_cds(mblocklist,graphid2ensemblid,ensemblid2graphid,cdsexon,cds2geneid,cds2geneexon)
    
    mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)
 
    if  compareExon == 'Yes':

        #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after completing cds\n")
        
        mblocklist =  merge_compatible_unordered(mblocklist,allcdsseq)
        
        #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after merge_compatible_unordered\n")

        mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)

        # i = 1
        # new_nb_mblock = len(mblocklist)
        # nb_mblock = new_nb_mblock + 1
        # while(new_nb_mblock != nb_mblock):
        #     print("in merge")
        mblocklist = merge_compatible_ordered(mblocklist,allcdsseq)
        #nb_mblock = new_nb_mblock
        #new_nb_mblock = len(mblocklist)
        #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after merge_compatible_ordered\n")
            # i+=1
        
        mblocklist =  merge_compatible_extremity(mblocklist)

        #print(time.time()-temps, "\n" ,len(mblocklist), "mblocks after merge_compatible_extremity\n")

        mblocklist =  order_mblocklist(mblocklist,graphid2ensemblid,ensemblid2graphid,gene2cds,cds2geneid)

    mblocklist = changesequenceid(mblocklist,graphid2ensemblid)

    print(time.time()-temps, "\n" ,len(mblocklist), "final mblocks")

    return mblocklist


