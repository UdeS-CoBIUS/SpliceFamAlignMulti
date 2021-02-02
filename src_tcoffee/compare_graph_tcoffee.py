#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``compare.py`` **module description**:

This module is the module that launches the comparison for all pairs of source CDS and target gene.

.. moduleauthor:: Aida Ouangraoua

2019 

"""

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import time
import os


def sort_block(targetdata, nbinitialsource, extendedsourcedata, comparisonresults,comparisonresults_idty,interval,cdsexon,cds2geneexon):
    graphid2ensemblid = {}   # dictionnary graphid to ensemblid
    ensemblid2graphid = {}   # dictionnary graphid to ensemblid
    allcdsseq = {}
    allgeneseq = {}
    identity_blocks = []
    interval_blocks = []
    sorted_blocks = []
    
    i = 0
    for gene in targetdata:
        geneid, geneseq = gene
        graphid2ensemblid["g"+str(i)] = geneid
        ensemblid2graphid[geneid] = "g"+str(i)
        allgeneseq[geneid] = geneseq
        i += 1

    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        graphid2ensemblid["c"+str(j)] = cdsid
        ensemblid2graphid[cdsid] = "c"+str(j)
        allcdsseq["c"+str(j)] = cdsseq

    i = 0
    for gene in targetdata:
        geneid,geneseq = gene
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            null,blocklist,null,null,null = comparisonresults[i][j]
            for k in range(len(blocklist)):
                block = blocklist[k]
                gstart,gend = block[2:]
                cstart, cend = block[:2]
                idty = comparisonresults_idty[i][j][k]
                length = cend-cstart
                if(cdsgeneid == geneid):
                    identity_blocks.append([geneid,cdsid,gstart,gend,cstart, cend, idty, length])
                else:
                    segment_cds = [cstart,cend]
                    segment_gene_cds = cds2genelocation(cdsid,segment_cds,cdsexon,cds2geneexon)
                    if(segment_gene_cds[1]-segment_gene_cds[0] == length):
                        sorted_blocks.append([geneid,cdsid,gstart,gend,cstart, cend, idty, length])
        i += 1
    sorted_blocks_interval = sort_block_interval(sorted_blocks,interval)
    sorted_blocks = identity_blocks + sorted_blocks_interval
    return sorted_blocks, graphid2ensemblid, ensemblid2graphid, allcdsseq, allgeneseq

def sort_block_interval(sorted_blocks,interval):
    sorted_blocks = sorted(sorted_blocks, key=lambda x:x[6],reverse = True)
    sorted_blocks_interval = []
    min = 1
    j = 0
    while(j < len(sorted_blocks)-1):
        max = min
        min = max - interval
        i = j
        while(j < len(sorted_blocks)-1 and sorted_blocks[j][6] >= min):
            j += 1
        sorted_blocks_interval += sorted(sorted_blocks[i:j], key=lambda x:x[7],reverse = True)
    return sorted_blocks_interval

def create_cc(blocks,ensemblid2graphid,graphid2ensemblid,gene2cds,cds2geneid,geneexon,cdsexon,cds2geneexon,geneexonstart,geneexonend):
    connectedcomponents = []
    supportcc = []
    nbgenecc = []
    completecc = []
    graphnode2cc = {}
    selectblocks = []
    for i in range(len(blocks)):
        block = blocks[i]
        geneid,cdsid,gstart,gend,cstart, cend, idty, length = block
        if(gstart not in geneexonstart[geneid] and len(geneexonstart[geneid]) > 0):
            closest_start = geneexonstart[geneid][0]
            dist = abs(closest_start-gstart)
            for s in geneexonstart[geneid]:
                if(abs(s-gstart) < dist):
                    closest_start = s
                    dist = abs(s-gstart)
            if(dist <= 10 and closest_start < gend):
                gstart = closest_start
        if(gend not in geneexonend[geneid] and len(geneexonend[geneid]) > 0):
            closest_end = geneexonend[geneid][0]
            dist = abs(closest_end-gend)
            for s in geneexonend[geneid]:
                if(abs(s-gend) < dist):
                    closest_end = s
                    dist = abs(s-gend)
            if(dist <= 10 and gstart < closest_end):
                gend = closest_end
            blocks[i][2]=gstart
            blocks[i][3]=gend
            
        if(([gstart,gend] in geneexon[geneid] or len(geneexon[geneid])==0 or [cstart,cend] in cdsexon[cdsid])):
            selectblocks.append(blocks[i])

    #print("nb selected blocks", len(selectblocks))
    for i in range(len(selectblocks)):
        block = selectblocks[i]
        geneid,cdsid,gstart,gend,cstart, cend, idty, length = block
        if(abs(gend-gstart-cend+cstart) > 50):
            continue
        genegraphid = ensemblid2graphid[geneid]
        cdsgraphid = ensemblid2graphid[cdsid]
        # create a node for gene segment
        genenode = genegraphid+"_"+str(gstart)+"_"+str(gend)
        genecc = -1
        if(genenode in graphnode2cc.keys() and graphnode2cc[genenode] not in completecc):
            genecc = graphnode2cc[genenode]
        # create a node for cds segment
        cdsnode = cdsgraphid+"_"+str(cstart)+"_"+str(cend)
        cdscc = -1
        segment_cds = [int(x) for x in cdsnode.split("_")[-2:]]
        segment_gene_cds = cds2genelocation(cdsid,segment_cds,cdsexon,cds2geneexon)
        gcstart,gcend = segment_gene_cds
        gene_cds_graphid = ensemblid2graphid[cds2geneid[cdsid]]
        gene_cds_node = gene_cds_graphid+"_"+str(gcstart)+"_"+str(gcend)
        if(cdsnode in graphnode2cc.keys() and graphnode2cc[cdsnode] not in completecc):
            cdscc = graphnode2cc[cdsnode]
        elif(gene_cds_node in graphnode2cc.keys() and graphnode2cc[gene_cds_node] not in completecc):
            cdscc = graphnode2cc[gene_cds_node]

        current = -1
        # if both segment are not already nodes of graph
        if(genecc == -1 and cdscc == -1):
            # create new connected component
            connectedcomponents.append({})
            connectedcomponents[-1][genenode] = [gstart,gend] 
            connectedcomponents[-1][cdsnode] = [cstart,cend] 
            connectedcomponents[-1][gene_cds_node] = [gcstart,gcend]
            graphnode2cc[genenode] = len(connectedcomponents)-1
            graphnode2cc[cdsnode] = len(connectedcomponents)-1
            graphnode2cc[gene_cds_node] = len(connectedcomponents)-1
            supportcc.append([block])
            nbgenecc.append(2)
            current = len(connectedcomponents)-1
            #print(1,len(connectedcomponents)-1,len(connectedcomponents[-1].keys()))
        # if cds segment is already a node
        elif(genecc == -1):
            connectedcomponents[cdscc][genenode] = [gstart,gend]
            connectedcomponents[cdscc][cdsnode] = [cstart,cend]
            graphnode2cc[genenode] = cdscc
            graphnode2cc[cdsnode] = cdscc
            supportcc[cdscc].append(block)
            nbgenecc[cdscc] += 1
            current = cdscc
            #print(2,cdscc,len(connectedcomponents[cdscc].keys()))

        # if gene segment is already a node 
        elif(cdscc == -1):
            connectedcomponents[genecc][cdsnode] = [cstart,cend]
            connectedcomponents[genecc][gene_cds_node] = [gcstart,gcend]
            graphnode2cc[cdsnode] = genecc
            graphnode2cc[gene_cds_node] = genecc
            supportcc[genecc].append(block)
            nbgenecc[genecc] += 1
            current = genecc
            #print(3,genecc,len(connectedcomponents[genecc].keys()))

        # if both segments are already nodes of the graph in different cc
        elif(cdscc != genecc):
            connectedcomponents[genecc][cdsnode] = [cstart,cend]
            for cnode in connectedcomponents[cdscc].keys():
                connectedcomponents[genecc][cnode]=connectedcomponents[cdscc][cnode]
                graphnode2cc[cnode] = genecc
            supportcc[genecc].append(block)
            supportcc[genecc] += supportcc[cdscc]
            nbgenecc[genecc] += nbgenecc[cdscc] 
            current = genecc
            connectedcomponents[cdscc] = {}
            supportcc[cdscc] = []
            #print(4,genecc,len(connectedcomponents[genecc].keys()))

        else:
            connectedcomponents[genecc][cdsnode] = [cstart,cend]
            graphnode2cc[cdsnode] = genecc
            supportcc[genecc].append(block)
            current = genecc
            #print(5,genecc,len(connectedcomponents[genecc].keys()))

        #if(nbgenecc[current] > len(geneexon.keys())):
        if(nbgenecc[current] > 50):
            completecc.append(current)

    # remove empty cc
    nb_cc = len(connectedcomponents)
    for i in range(nb_cc - 1,-1,-1):
        if(len(connectedcomponents[i].keys())==0):
            del(connectedcomponents[i])
            del(supportcc[i])
        
    return connectedcomponents,supportcc

def generate_segment_matches_from_segments(segments, targetdata, genesegment,cds2geneid,cdsexon,cds2geneexon,idty_threshold):
    segment_matches = []
    for i in range(len(targetdata)):
        segment_matches.append([])
        for j in range(len(targetdata)):
            segment_matches[i].append([])

    gene_index = {}
    for i in range(len(targetdata)):
        gene_index[targetdata[i][0]] = i

    for segment in segments:
        geneid, cdsid, startgene, startcds, length, sim = segment
        genesegment = add_genesegment(genesegment, geneid,[startgene,startgene+length])
        cdsgeneid = cds2geneid[cdsid]
        segment_cds = [startcds,startcds+length]
        segmentgene_cds = cds2genelocation(cdsid,segment_cds,cdsexon,cds2geneexon)
        if(segment_cds[1]-segment_cds[0] == segmentgene_cds[1]-segmentgene_cds[0] and sim >= idty_threshold):
            if(gene_index[geneid] < gene_index[cdsgeneid]):
                segment_matches[gene_index[geneid]][gene_index[cdsgeneid]].append([startgene, segmentgene_cds[0],length, 1, sim*1000])
            elif(gene_index[cdsgeneid] < gene_index[geneid]):
                segment_matches[gene_index[cdsgeneid]][gene_index[geneid]].append([segmentgene_cds[0],startgene,length, 1, sim*1000])
    return segment_matches, genesegment

def generate_segment_matches(connectedcomponents, supportcc, targetdata, genesegment,graphid2ensemblid,ensemblid2graphid,cds2geneid,cdsexon,cds2geneexon,blockalignment):
    segment_matches = []
    for i in range(len(targetdata)):
        segment_matches.append([])
        for j in range(len(targetdata)):
            segment_matches[i].append([])

    gene_index = {}
    gene_sequence = {}
    for i in range(len(targetdata)):
        gene_index[targetdata[i][0]] = i
        gene_sequence[targetdata[i][0]] = targetdata[i][1]

    geneids = []
    for i in range(len(connectedcomponents)):
        geneids.append([])
        cci = connectedcomponents[i]
        for id in cci.keys():
            if id[0] == 'g':
                geneids[i].append(id)
                
    discard = set()
    for i in range(len(connectedcomponents)):
        cci = connectedcomponents[i]
        if(len(geneids[i]) <= 1):
            discard.add(i)
    #print("nb discarded cc", len(discard))

    for i in range(len(connectedcomponents)):
        if(i in discard):
            continue
        #print(i)
        seq_len = []
        cc = connectedcomponents[i]
        #print(i,len(cc.keys()))
        seqfile = "input.fasta"
        outputntfile = "outputnt.fasta"
        outputaafile = "outputaa.fasta"
        file = open(seqfile,"w")
        for id in geneids[i]:
            s1,e1 = id.split("_")[-2:]
            gid = id.split("_"+s1+"_"+e1)[0]
            genesegment = add_genesegment(genesegment, graphid2ensemblid[gid],cc[id])
            s1,e1 = cc[id]
            seq_len.append(e1-s1)
            file.write(">"+graphid2ensemblid[gid]+"_"+str(s1)+"_"+str(e1) + "\n" + gene_sequence[graphid2ensemblid[gid]][s1:e1]+"\n")
        file.close()
        align = ""
        if(min(seq_len) == max(seq_len)):
            align = AlignIO.read(seqfile, "fasta")

        else:
            if(blockalignment == "mafft"):
                mafft_cline = MafftCommandline(input=seqfile,localpair=True,lop = -6.0)
                stdout, stderr = mafft_cline()
                align = AlignIO.read(StringIO(stdout), "fasta")
            elif(blockalignment == "macse"):
                os.system('java -jar src_multi/macse_v2.03.jar -prog alignSequences -seq '+  seqfile +' -out_NT '+ outputntfile +' -out_AA '+outputaafile+" 2>/dev/null"+" >/dev/null")
                align = AlignIO.read(outputntfile, "fasta")
                
        alignment = {}
        alnsegment = {}
        for record in align:
            #print(record.id)
            alignment[record.id] = record.seq
            alnsegment[record.id] = []
            k = 0
            s = -1
            e = -1
            while(k < len(record.seq)):
                if(s == -1 and record.seq[k] != '-' and record.seq[k] != '!'):
                    s = k
                elif(s != -1 and e == -1 and record.seq[k] == '-' and record.seq[k] != '!'):
                    e = k
                    alnsegment[record.id].append([s,e])
                    s = -1
                    e = -1
                k += 1
            if(s != -1 and e == -1):
               e = k
               alnsegment[record.id].append([s,e])
        #print()       
        support = supportcc[i]
        genesegmentpairs = {}
        for block in support:
            geneid,cdsid,gstart,gend,cstart, cend, idty, length = block
            geneseq = geneid+"_"+str(gstart)+"_"+str(gend)
            cdsseq = cdsid+"_"+str(cstart)+"_"+str(cend)
            segment_cds = [int(x) for x in cdsseq.split("_")[-2:]]
            segment_gene_cds = cds2genelocation(cdsid,segment_cds,cdsexon,cds2geneexon)
            gcstart,gcend = segment_gene_cds
            gene_cds_seq = cds2geneid[cdsid]+"_"+str(gcstart)+"_"+str(gcend)

            pairid = geneseq+":"+gene_cds_seq
            #if(pairid in genesegmentpairs.keys()):
            genesegmentpairs[pairid] = int(idty*1000)
            #else:
                #genesegmentpairs[pairid] = 0
        for pairid in genesegmentpairs.keys():
            id1,id2 = pairid.split(":")
            #print(id1,id2)
            union = unionsegment(sorted(alnsegment[id1]+alnsegment[id2]))
            unionlen = sum([e-s for s,e in union])
            s1,e1 = id1.split("_")[-2:]
            gid1 = id1.split("_"+s1+"_"+e1)[0]
            s2,e2 = id2.split("_")[-2:]
            gid2 = id2.split("_"+s2+"_"+e2)[0]
            idty = genesegmentpairs[pairid]
            if(gene_index[gid1] < gene_index[gid2]):
                common = commonsegment(alnsegment[id1],alnsegment[id2])
                commonlen = sum([e-s for s,e,null,null in common])
                sim = int(1000.0*commonlen/unionlen)
                if(sim > 700):
                    for s,e,l1,l2 in common:
                        segment_matches[gene_index[gid1]][gene_index[gid2]].append([cc[ensemblid2graphid[gid1]+"_"+str(s1)+"_"+str(e1)][0]+l1, cc[ensemblid2graphid[gid2]+"_"+str(s2)+"_"+str(e2)][0]+l2,e-s, 1, idty])
            else:
                common = commonsegment(alnsegment[id2],alnsegment[id1])
                commonlen = sum([e-s for s,e,null,null in common])
                sim = int(1000.0*commonlen/unionlen)
                if(sim > 700):
                    for s,e,l2,l1 in common:
                        segment_matches[gene_index[gid2]][gene_index[gid1]].append([cc[ensemblid2graphid[gid2]+"_"+str(s2)+"_"+str(e2)][0]+l2, cc[ensemblid2graphid[gid1]+"_"+str(s1)+"_"+str(e1)][0]+l1,e-s, 1, idty])
    return segment_matches, genesegment

def write_segment_matches(targetdata, nbinitialsource, extendedsourcedata, segment_matches, genesegment,cdsexon, cds2geneexon, fasta_filename, output_filename,outputalignment):
    genelength = []
    geneindex = {}
    output_fasta = open(fasta_filename,"w")
    output_segment = open(output_filename,"w")
    if(outputalignment == "tcoffee"):
        output_segment.write("! TC_LIB_FORMAT_01\n")
        output_segment.write(str(len(targetdata)+nbinitialsource)+"\n")
    for i in range(len(targetdata)):
        geneid,geneseq = targetdata[i]
        geneindex[geneid] = i
        genesegment[geneid] = sorted(genesegment[geneid])
        geneseq_ = ""
        for segment in genesegment[geneid]:
            geneseq_ += geneseq[segment[0]:segment[1]]
        if(outputalignment == "tcoffee"):
            output_segment.write(geneid + " " + str(len(geneseq_)) + " " + geneseq_ + "\n")
        output_fasta.write(">" + geneid + "\n" + geneseq_ + "\n")
        genelength.append(len(geneseq_))
    for j in range(nbinitialsource):
        cds = extendedsourcedata[j]
        cdsid,cdsseq,cdsgeneid,null = cds
        if(outputalignment == "tcoffee"):
            output_segment.write(cdsid + " " + str(len(cdsseq)) + " " + cdsseq + "\n")
        output_fasta.write(">" + cdsid + "\n" + cdsseq + "\n")
    output_fasta.close()

    segment_matches_ = []
    for i in range(len(targetdata)):
        segment_matches_.append([])
        for j in range(len(targetdata)):
            segment_matches_[i].append([])

    for i in range(len(targetdata)):
        for j in range(i+1,len(targetdata)):
            for sm in segment_matches[i][j]:
                start1,start2,length,lsim,gsim = sm
                overlap = []
                for sm_ in segment_matches_[i][j]:
                    start1_,start2_,length_,lsim_,gsim_ = sm_
                    if(start1_ - start1 == start2_ - start2 and start1 + length > start1_ and start1_ + length_ > start1):
                        overlap.append(sm_)
                if(len(overlap) > 0):
                    end1 = max([(x[0]+x[2]) for x in overlap]+[start1+length])
                    start1 = min([x[0] for x in overlap]+[start1])
                    start2 = min([x[1] for x in overlap]+[start2])
                    length = end1 - start1
                    gsim = min([x[4] for x in overlap]+[gsim])
                segment_matches_[i][j].append([start1,start2,length,lsim,gsim])
                for sm_ in overlap:
                    segment_matches_[i][j].remove(sm_)
                    
    for i in range(len(targetdata)):
        geneid = targetdata[i][0]
        if(outputalignment == "seqan_tcoffee"):
            output_segment.write(">"+geneid+"\n")
        for j in range(i+1,len(targetdata)):
            if(outputalignment == "tcoffee"):
                output_segment.write("#" + str(i+1) + " " + str(j+1)+"\n")
                
            for sm in segment_matches_[i][j]:
                start1,start2,length,lsim,gsim = sm
                start1_ = 0
                k = 0
                while(start1 >= genesegment[targetdata[i][0]][k][1]):
                    start1_ += genesegment[targetdata[i][0]][k][1] - genesegment[targetdata[i][0]][k][0]
                    k += 1
                start1_ += start1 - genesegment[targetdata[i][0]][k][0]
                start2_ = 0
                k = 0
                while(start2 >= genesegment[targetdata[j][0]][k][1]):
                    start2_ += genesegment[targetdata[j][0]][k][1] - genesegment[targetdata[j][0]][k][0]
                    k += 1
                start2_ += start2 - genesegment[targetdata[j][0]][k][0]
                #assert(start1_+length<=genelength[i])
                #assert(start2_+length<=genelength[j])
                
                if(outputalignment == "tcoffee"):
                    output_segment.write("+BLOCK+ " + str(length) + " " + str(start1_+1) + " " + str(start2_+1) + " " + str(gsim)+"\n")
                elif(outputalignment == "seqan_tcoffee"):
                    output_segment.write(" "+targetdata[j][0]+"\t"+ str(start1_+1) + "\t" + str(start2_+1) +  "\t" +str(length)+"\n")
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            if(cdsgeneid != geneid):
                continue
            if(outputalignment == "tcoffee"):
                output_segment.write("#" + str(i+1) + " " + str(len(targetdata)+j+1)+"\n")
        
            for exon in cdsexon[cdsid]:
                genestart = cds2geneexon[cdsid][exon[0]][0]
                genestart_ = 0
                k = 0
                while(genestart >= genesegment[cdsgeneid][k][1]):
                    genestart_ += genesegment[cdsgeneid][k][1] - genesegment[cdsgeneid][k][0]
                    k += 1
                genestart_ += genestart - genesegment[cdsgeneid][k][0]

                if(outputalignment == "tcoffee"):
                    output_segment.write("+BLOCK+ " + str(exon[1]-exon[0]) + " " + str(genestart_+1) + " " + str(exon[0]+1) + " " + str(10000)+"\n")
                elif(outputalignment == "seqan_tcoffee"):
                    output_segment.write(" "+cdsid+"\t"+ str(genestart_+1) + "\t" + str(exon[0]+1) +  "\t" +str(exon[1]-exon[0])+"\n")
    if(outputalignment == "tcoffee"):
        output_segment.write("! SEQ_1_TO_N\n")                
    output_segment.close()

def unionsegment(alnsegment):
    union = []
    if(len(alnsegment) > 0):
        s,e = alnsegment[0]
        union.append([s,e])
        for i in range(1, len(alnsegment)):
            s,e = alnsegment[i]
            if(s <= union[-1][1]):
                union[-1][1] = max(union[-1][1],e)
            else:
                union.append([s,e])
    return union
       
def commonsegment(alnsegment1,alnsegment2):
    common = []
    i1 = 0
    i2 = 0
    len1 = 0
    len2 = 0
    while(i1 < len(alnsegment1) and i2 < len(alnsegment2)):
        s1,e1 = alnsegment1[i1]
        s2,e2 = alnsegment2[i2]
        if(s2 < s1):
            if(s1 < e2):
                if(e1 <= e2):
                    common.append([s1,e1, len1, len2 + s1-s2])
                    i1 += 1
                    len1 += e1-s1
                else:
                    common.append([s1,e2, len1, len2 + s1-s2])
                    i2 += 1
                    len2 += e2-s2
            else:
                i2 += 1
                len2 += e2-s2
        elif(s1 < s2):
            if(s2 < e1):
                if(e2 <= e1):
                    common.append([s2,e2, len1+s2-s1, len2])
                    i2 += 1
                    len2 += e2-s2
                else:
                    common.append([s2,e1, len1+s2-s1, len2])
                    i1 += 1
                    len1 += e1-s1
            else:
                i1 += 1
                len1 += e1-s1
        else:
            if(e1 <= e2):
                common.append([s1,e1, len1, len2])
                i1 += 1
                len1 += e1-s1
            else:
                common.append([s2,e2, len1, len2])
                i2 += 1
                len2 += e2-s2
    return common
    
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


def order_mblocklist(mblocklist):
    """
    This function orders blocks in the blocklist by increasing
    order of query start location.
    """
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


def add_genesegment(genesegment, geneid, exon):
    if(geneid in (genesegment.keys())):
        overlap = []
        for i in range(len(genesegment[geneid])):
            segment = genesegment[geneid][i]
            if(segment[1] > exon[0] and exon[1] > segment[0]):
                overlap.append(segment)
        if(len(overlap) == 0):
            genesegment[geneid].append(exon)
        elif(len(overlap) == 1):
            segment1 = [min(overlap[0][0],exon[0]), overlap[0][0]]
            segment2 = [overlap[0][1], max(overlap[0][1],exon[1])]
            if(segment1[0] < segment1[1]):
                genesegment[geneid].append(segment1)
            if(segment2[0] < segment2[1]):
                genesegment[geneid].append(segment2)
        else:
            segment1 = [min([x[0] for x in overlap]+[exon[0]]), min([x[0] for x in overlap])]
            segment2 = [max([x[1] for x in overlap]), max([x[1] for x in overlap]+[exon[1]])]
            if(segment1[0] < segment1[1]):
                genesegment[geneid].append(segment1)
            if(segment2[0] < segment2[1]):
                genesegment[geneid].append(segment2)
            overlap=sorted(overlap)
            for j in range(len(overlap)-1):
                segment = [overlap[j][1],overlap[j+1][0]]
                if(segment[0] < segment[1]):
                    genesegment[geneid].append(segment)
    else:
        genesegment[geneid] = [exon]
    genesegment[geneid] = sorted(genesegment[geneid])
    return genesegment

def compute_geneexon(geneexon):
    geneexonstart = {}
    geneexonend = {}
    genesegment = {}
    for geneid in geneexon.keys():
        geneexonstart[geneid] = sorted(list(set([x[0] for x in geneexon[geneid]])))
        geneexonend[geneid] = sorted(list(set([x[1] for x in geneexon[geneid]])))
        genesegment[geneid] = []
        for exon in geneexon[geneid]:
            overlap = []
            for i in range(len(genesegment[geneid])):
                segment = genesegment[geneid][i]
                if(segment[1] > exon[0] and exon[1] > segment[0]):
                    overlap.append(segment)
            if(len(overlap) == 0):
                genesegment[geneid].append(exon)
            else:
                segment = [min([x[0] for x in overlap]+[exon[0]]), max([x[1] for x in overlap]+[exon[1]])]
                for s in overlap:
                    genesegment[geneid].remove(s)     
                genesegment[geneid].append(segment)
        genesegment[geneid] = sorted(genesegment[geneid])
    return geneexonstart, geneexonend, genesegment

def compute_msa(extendedsourcedata,targetdata,comparisonresults,comparisonresults_idty,geneexon,cdsexon,nbinitialsource,cds2geneid,cds2geneexon,gene2cds, segments, idty_threshold, compareExon, outputprefix,blockalignment,outputalignment):

    temps=time.time()

    geneexonstart, geneexonend, genesegment = compute_geneexon(geneexon)

    genesegment_init = genesegment.copy()

    segment_matches = []
    if(len(segments) == 0):
        interval = 0.01 
        blocks, graphid2ensemblid, ensemblid2graphid, allcdsseq, allgeneseq = sort_block(targetdata, nbinitialsource, extendedsourcedata, comparisonresults, comparisonresults_idty,interval, cdsexon,cds2geneexon)
        #print("nb blocks: ", len(blocks))
        connectedcomponents,supportcc = create_cc(blocks,ensemblid2graphid,graphid2ensemblid,gene2cds,cds2geneid,geneexon,cdsexon,cds2geneexon,geneexonstart,geneexonend)
        #print("initial nb cc: ", len(connectedcomponents))

        segment_matches,genesegment = generate_segment_matches(connectedcomponents, supportcc, targetdata, genesegment,graphid2ensemblid,ensemblid2graphid,cds2geneid,cdsexon,cds2geneexon,blockalignment)
    else:
        segment_matches,genesegment =  generate_segment_matches_from_segments(segments, targetdata, genesegment,cds2geneid,cdsexon,cds2geneexon,idty_threshold)
    #print(time.time()-temps, "segment_matches")

    temps=time.time()
    tcoffee(targetdata, nbinitialsource, extendedsourcedata, segment_matches, genesegment,genesegment_init,cdsexon, cds2geneexon, outputprefix,outputalignment)
    #print(time.time()-temps, "tcofeee")

    mblocklist= []
    #mblocklist= initialize_mblocklist(connectedcomponents, graphid2ensemblid)
    
    #mblocklist = changesequenceid(mblocklist,graphid2ensemblid)

    #print("before order")
    mblocklist =  order_mblocklist(mblocklist)
    #print(time.time()-temps, "\n" ,len(mblocklist), "final mblocks")

    
    return mblocklist

def tcoffee(targetdata, nbinitialsource, extendedsourcedata, segment_matches, genesegment,genesegment_init,cdsexon, cds2geneexon, outputprefix,outputalignment):
    fasta_filename = outputprefix+"genesegment.fasta"
    segment_filename = outputprefix+"segment_matches"
    if(outputalignment == "tcoffee"):
        segment_filename += ".tc_lib"
    elif(outputalignment == "seqan_tcoffee"):
        segment_filename += ".mums"
    aln_filename1 = outputprefix+"genesegmentaln.aln"
    aln_filename2 = outputprefix+"genesegmentaln.fasta"
        
    write_segment_matches(targetdata, nbinitialsource, extendedsourcedata, segment_matches, genesegment,cdsexon, cds2geneexon, fasta_filename, segment_filename,outputalignment)

    align = ""
    if(outputalignment == "tcoffee"):
        t_coffee_command = "t_coffee -in " + fasta_filename +  " -lib " +  segment_filename +  " -outfile "+aln_filename1+" 2>/dev/null"+" >/dev/null"
        os.system(t_coffee_command)
        #t_coffee_command = "t_coffee -other_pg seq_reformat -in " + aln_filename1 + " -output fasta_aln > " + aln_filename2 
        #os.system(t_coffee_command)
        align = AlignIO.read(aln_filename1, "clustal")
    elif(outputalignment == "seqan_tcoffee"):
        seqan_coffee_command = "seqan_tcoffee -s " + fasta_filename +  " -l " +  segment_filename +  " -o "+aln_filename2+" 2>/dev/null"+" >/dev/null"
        os.system(seqan_coffee_command)
        
        align = AlignIO.read(aln_filename2, "fasta")
    alignment = {}
    length = 0
    for record in align:
        alignment[record.id] = record.seq
        length = len(record.seq)

    #print(length)
    removecolumn = []
    nbremoved = 0
    for i in range(length):
        if(all([alignment[id][i] == '-' for id in cdsexon.keys()])):
            nbremoved += 1
            if(len(removecolumn)== 0):
                removecolumn.append([i,i+1])
            elif(i == removecolumn[-1][1]):
                removecolumn[-1][1] += 1
            else:
                removecolumn.append([i,i+1])
    #print(nbremoved)
    #print(removecolumn)
    
    # nbredisue = {}
    # for id in alignment.keys():
    #     nbredisue[id] = 0
    # for i in range(length):
    #     for id in alignment.keys():
    #         if(alignment[id][i] != '-'):
    #             nbredisue[id] += 1
    #         if(id in cdsexon.keys() and (nbredisue[id]-1) in )


