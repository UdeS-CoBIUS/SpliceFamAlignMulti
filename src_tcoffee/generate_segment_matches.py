#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``generate_segment_matches.py`` **module description**:

This module generate a set of segment matches given a pairwise spliced alignments of CDS and genes of a gene family

.. moduleauthor:: Aida Ouangraoua

2020

"""

from Bio import SeqIO
from Bio import AlignIO
from Bio import pairwise2
import argparse
import time
import os


def get_data_from_files(args):
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
    comparisonresults_length		=	[]
    geneexon = {}
    cdsexon = {}
    liste_CDS_GENE		=	[]
    liste_CDS_GENE_idty		=	[]
    liste_CDS_GENE_length		=	[]
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
            CDS_GENE_length = []
            while i< len(lines) and not(lines[i] in ['\n', '\r\n']):
                line = lines[i]	
                tab = line.split("\t")
                if len(tab) == 9:
                    idCDS			= 	tab[0]
                    idGene			= 	tab[1]
                    length		= 	int(tab[2])
                    beginCDS		= 	int(tab[3])
                    endCDS			= 	int(tab[4])
                    beginGene		=	int(tab[5])
                    endGene			=	int(tab[6])
                    part_CDS_Gene 	= 	[beginCDS, endCDS, beginGene, endGene]
                    part_CDS_Gene_idty = float(tab[7])
                    part_CDS_Gene_length = int(tab[2])
                    # if(idGene == cds2gene[idCDS]):
                    #     if(part_CDS_Gene_idty != 1.0):
                    #         print(idCDS,idGene)

                    if(part_CDS_Gene_idty >= idty_threshold and 0 <= beginCDS <= endCDS and 0 <= beginGene <= endGene):
                        CDS_GENE[1].append(part_CDS_Gene)
                        CDS_GENE_idty.append(part_CDS_Gene_idty)
                        CDS_GENE_length.append(part_CDS_Gene_length)
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
                liste_GENE_length.append(CDS_GENE_length)
            else:
                if(len(liste_GENE) != 0):
                    comparisonresults.append(liste_GENE)
                    comparisonresults_idty.append(liste_GENE_idty)
                    comparisonresults_length.append(liste_GENE_length)
                liste_GENE    =	[CDS_GENE]
                liste_GENE_idty    =	[CDS_GENE_idty]
                liste_GENE_length    =	[CDS_GENE_length]
            prev_GENE = idGene
        else:
            i += 1

    comparisonresults.append(liste_GENE)
    comparisonresults_idty.append(liste_GENE_idty)
    comparisonresults_length.append(liste_GENE_length)
    return comparisonresults,comparisonresults_idty,comparisonresults_length,geneexon,cdsexon,cds2gene,cds2geneexon,gene2cds

def generate_segment_matches(targetdata, nbinitialsource, extendedsourcedata, comparisonresults,comparisonresults_idty,cdsexon,cds2geneexon,cds2geneid,genesegment):
    segment_matches = []
    for i in range(len(targetdata)):
        segment_matches.append([])
        for j in range(len(targetdata)):
            segment_matches[i].append([])

    gene_index = {}
    for i in range(len(targetdata)):
        gene_index[targetdata[i][0]] = i
    
    segment_matches_id = []
    exon_matches_id = []
    ii = 0
    for gene in targetdata:
        geneid,geneseq = gene
        for j in range(nbinitialsource):
            cds = extendedsourcedata[j]
            cdsid,cdsseq,cdsgeneid,null = cds
            exonendlist = [x[1] for x in cdsexon[cdsid]]
            if(geneid !=  cds2geneid[cdsid]):
                null,blocklist,null,null,null = comparisonresults[ii][j]
                blocklist_idty = comparisonresults_idty[ii][j]
                for k in range(len(blocklist)):
                    block = blocklist[k]
                    idty = int(blocklist_idty[k]*1000)
                    gstart,gend = block[2:]
                    cstart, cend = block[:2]
                    gcstart, gcend = cds2genelocationinit(cdsid,[cstart,cend],cdsexon,cds2geneexon)
                    id = ""
                    if(gene_index[geneid] < gene_index[cds2geneid[cdsid]]):
                        id = geneid + "_" + cds2geneid[cdsid] + "_" + str(gstart)+ "_" + str(gend)+ "_" + str(gcstart)+ "_" + str(gcend)
                    else:
                        id = cds2geneid[cdsid] + "_" + geneid + "_" + str(gcstart)+ "_" + str(gcend)+ "_" + str(gstart_)+ "_" + str(gend)
        
                    if(id not in exon_matches_id):
                        exon_matches_id.append(id)

                        genesegment = add_genesegment(genesegment, geneid, [gstart,gend])
                        gseq = geneseq[gstart:gend]
                        cseq = cdsseq[cstart:cend]
                        gseq_=""
                        cseq_=""
                        if(gend-gstart == cend-cstart):
                            gseq_=gseq
                            cseq_=cseq
                        else:
                            # maximise matches and minimize gaps
                            alignment = pairwise2.align.globalms(gseq, cseq,2,0,-5,-1)
                            gseq_,cseq_= alignment[0][0],alignment[0][1]
                        length = len(cseq_)
                        i = 0
                        ig = 0
                        ic = 0
                        while(i < length):
                            while(i < length and (gseq_[i] == '-' or cseq_[i] == '-')):
                                if(gseq_[i] != '-'):
                                    ig += 1
                                if(cseq_[i] != '-'):
                                    ic += 1
                                i+=1
                            start_ = i
                            cstart_ = ic
                            gstart_ = ig
                            match = 0
                            while(i < length and gseq_[i] != '-' and cseq_[i] != '-' and ((i-start_ == 0) or (cstart+ic) not in exonendlist)):
                                ig += 1
                                ic += 1
                                if(cseq_[i] == gseq_[i]):
                                    match += 1
                                i+=1
                            end_ = i
                            length_ = end_-start_
                            if(length_ > 0):
                                sim = 1.0 *(match/length_)
                                gcstart_,gcend_ = cds2genelocation(cdsid,[cstart+cstart_,cstart+cstart_+length_],cdsexon,cds2geneexon)
                                id = ""
                                if(gene_index[geneid] < gene_index[cds2geneid[cdsid]]):
                                    id = geneid + "_" + cds2geneid[cdsid] + "_" + str(gstart+gstart_)+ "_" + str(gcstart_)+ "_" + str(length_)
                                else:
                                    id = cds2geneid[cdsid] + "_" + geneid + "_" + str(gcstart_)+ "_" + str(gstart+gstart_)+ "_" + str(length_)
                
                                if(id not in segment_matches_id):
                                    segment_matches_id.append(id)
                                    if(gene_index[geneid] < gene_index[cds2geneid[cdsid]]):
                                        segment_matches[gene_index[geneid]][gene_index[cds2geneid[cdsid]]].append([gstart+gstart_, gcstart_, length_, sim, idty])
                                    else:
                                        segment_matches[gene_index[cds2geneid[cdsid]]][gene_index[geneid]].append([gcstart_, gstart+gstart_, length_, sim, idty])

        ii += 1
    return segment_matches,genesegment

def write_segment_matches(targetdata, segment_matches, genesegment, fasta_filename, output_filename):
    genelength = []
    output_fasta = open(fasta_filename,"w")
    output_segment = open(output_filename,"w")
    output_segment.write("! TC_LIB_FORMAT_01\n")
    output_segment.write(str(len(targetdata))+"\n")
    for i in range(len(targetdata)):
        geneid,geneseq = targetdata[i]
        genesegment[geneid] = sorted(genesegment[geneid])
        geneseq_ = ""
        for segment in genesegment[geneid]:
            geneseq_ += geneseq[segment[0]:segment[1]]
        output_segment.write(geneid + " " + str(len(geneseq_)) + " " + geneseq_ + "\n")
        output_fasta.write(">" + geneid + "\n" + geneseq_ + "\n")
        genelength.append(len(geneseq_))
    output_fasta.close()
    for i in range(len(targetdata)):
        for j in range(i+1,len(targetdata)):
            output_segment.write("#" + str(i+1) + " " + str(j+1)+"\n")
            for sm in segment_matches[i][j]:
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
                assert(start1_+length<=genelength[i])
                assert(start2_+length<=genelength[j])
                
                output_segment.write("+BLOCK+ " + str(length) + " " + str(start1_+1) + " " + str(start2_+1) + " " + str(gsim)+"\n")
    output_segment.write("! SEQ_1_TO_N\n")                
    output_segment.close()


                          
# Compute the gene location given the cds location
def cds2genelocation(cdsid,segment,cdsexon,cds2geneexon):
    start,end = segment
    gstart,gend = 0,0
    exonstart = -1
    exonend = -1
    loc = []

    for i in range(len(cdsexon[cdsid])):
        exon = cdsexon[cdsid][i]
        if(exon[0] <= start and start < exon[1]):
            gstart = cds2geneexon[cdsid][exon[0]][0] + (start-exon[0])
            assert(exon[0] < end and end <= exon[1])
            gend = cds2geneexon[cdsid][exon[0]][0] + (end-exon[0])

    return [gstart,gend]

def cds2genelocationinit(cdsid,segment,cdsexon,cds2geneexon):
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

def compute_genesegment(geneexon):
    genesegment = {}
    for geneid in geneexon.keys():
        genesegment[geneid] = []
        for exon in geneexon[geneid]:
            add_genesegment(genesegment, geneid, exon)
        genesegment[geneid] = sorted(genesegment[geneid])
    return genesegment

def add_genesegment(genesegment, geneid, exon):
    if(geneid in (genesegment.keys())):
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
    else:
        genesegment[geneid] = [exon]
    genesegment[geneid] = sorted(genesegment[geneid])
    return genesegment

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

    parser.add_argument('-ce', '--compareExon', help="The method includes a final step that compares exons of blocks for further merges: Yes or No")

    parser.add_argument('-msa', '--msaMethod', help="Multiple sequence aligner: muscle or mafft")
    
    return parser

def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    treefile = args.treeFile
    idty_threshold = float(args.identityThreshold)
    outputprefix = args.outputPrefix
    pairwisealnfile = args.pairwiseAlnFile

    if (outputprefix == None):
        print("Argument -op <outputPrefix> is required")
        return
            
    if (pairwisealnfile == None):
        print("Argument -palnf <pairwiseAlnFile> is required")
        return

    if(outputprefix != None and pairwisealnfile != None):
        print("Retrieving input data...")
        sourcedata,targetdata = get_data_from_files(args)
        nbinitialsource = len(sourcedata)

        print("Reading pairwise alignment file...")
        comparisonresults,comparisonresults_idty,comparisonresults_length,geneexon,cdsexon,cds2geneid,cds2geneexon,gene2cds = parseResultFile(pairwisealnfile,sourcedata, idty_threshold)
        
        temps=time.time()
        print("Computing multiple alignement...")

        genesegment = compute_genesegment(geneexon)
        genesegment_init = genesegment.copy()
        print(time.time()-temps)

        segment_matches,genesegment = generate_segment_matches(targetdata, nbinitialsource, sourcedata, comparisonresults,comparisonresults_idty,cdsexon,cds2geneexon,cds2geneid,genesegment)
        print(time.time()-temps)

        segment_filename = outputprefix+"segment_matches.tc_lib"
        fasta_filename = outputprefix+"genesegment.fasta"
        write_segment_matches(targetdata, segment_matches, genesegment, fasta_filename, segment_filename)
        print(time.time()-temps)

        t_coffee_command = "t_coffee -in " + fasta_filename +  " -lib " +  segment_filename +  " -outfile test.aln"
        os.system(t_coffee_command)
        print(time.time()-temps)

if __name__ == '__main__':
    main()
