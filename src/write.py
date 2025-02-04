#!/usr/bin/python2.7
#-*- coding: utf-8 -*-

"""

``write.py`` **module description**:

This module is the module that formats and writes the results in output files.

.. moduleauthor:: Safa Jammali 

2017-2018

"""

import multiprocessing
from functools import partial
from contextlib import contextmanager
from multiprocessing import Pool
from utils import *
import pathlib
sys.path.append(str(os.path.abspath(os.path.dirname(sys.argv[0])))+ '/fsepsa')
import scoring_matrix
from scoring_matrix import *
#import fse
from fse import fse
import translator
from translator import aamap

BORDER_LENGTH = 80
LENGTH_ALN_LINE = 80

#############################
##  INPUT FILE WRITTING #####
############################

def  write_input_files(outputprefix,sourcedata,targetdata):

    initsourcefile = open(outputprefix+"initialsource.fasta","w")
    targetfile = open(outputprefix+"target.fasta","w")
    source2targetfile = open(outputprefix+"initialsource2target.txt","w")
    initsourceexonlistfile = open(outputprefix+"initialsourceexonlist.txt","w")

    for cds in sourcedata:
        cdsid,cdsseq,cdsgeneid,exonlist = cds
	
        initsourcefile.write(">"+cdsid+"\n")
        initsourcefile.write(cdsseq+"\n\n")
        source2targetfile.write(cdsid + " " + cdsgeneid + "\n")
        initsourceexonlistfile.write(">"+cdsid+"\n")
        for block in exonlist:
            cdsstart,cdsend,genestart,geneend = block
            initsourceexonlistfile.write(str(cdsstart) + " " + str(cdsend) + " " + str(genestart) + " " + str(geneend)+"\n")
        initsourceexonlistfile.write("\n")
    initsourcefile.close()
    source2targetfile.close()
    initsourceexonlistfile.close()
    for gene in targetdata:
        geneid,geneseq = gene
        targetfile.write(">"+geneid+"\n")
        targetfile.write(geneseq+"\n\n")
    targetfile.close()

#############################
## OUTPUT FILE WRITTING #####
############################
        
def writeOutfile(outputprefix,outputformat,outputalignment,extendedsourcedata,targetdata,extendedorthologygroups,comparisonresults,nbinitialsource):
    

    extendedsourcefile_buffer = ""
    extendedsource2targetfile_buffer = ""
    for cds in extendedsourcedata:
        cdsid,cdsseq,cdsgeneid ,null = cds 
        extendedsourcefile_buffer += ">"+cdsid+"\n"
        extendedsourcefile_buffer += cdsseq+"\n\n"
        extendedsource2targetfile_buffer += cdsid + " " + cdsgeneid + "\n"

    extendedsourcefile = open(outputprefix+"source.fasta","w")
    extendedsourcefile.write(extendedsourcefile_buffer)
    extendedsourcefile.close()

    extendedsource2targetfile = open(outputprefix+"source2target.txt","w")
    extendedsource2targetfile.write(extendedsource2targetfile_buffer)
    extendedsource2targetfile.close()

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
    pool_results = p.map(partial(pool_print_blocklist, outputformat = outputformat, outputalignment = outputalignment), datalist)

    #pool_results = []
    #for x in datalist:
        #pool_results.append(pool_print_blocklist(x,outputformat))
        
    results = {}
    for item in pool_results:
        results[item[0] + item[1]] = [item[2], item[3], item[4]]
            
    resultfile_buffer = ""
    cdsexonendlistfile_buffer = ""
    segmentfile_buffer = ""
    
    for gene in targetdata:
        geneid,geneseq = gene
        for i in range(nbinitialsource):
            cds = extendedsourcedata[i]
            cdsid,cdsseq,cdsgeneid ,null = cds
            resultfile_buffer += results[cdsid + geneid][0]+"\n"
            cdsexonendlistfile_buffer += results[cdsid + geneid][1]+"\n"
            for startgene, startcds, length, sim in results[cdsid + geneid][2]:
                segmentfile_buffer += geneid + "\t" + cdsid + "\t" + str(startgene)+ "\t" + str(startcds)+ "\t" + str(length)+ "\t" + sim+"\n"
            
    resultfile = open(outputprefix+"result.txt","w")
    resultfile.write(resultfile_buffer )
    resultfile.close()

    cdsexonendlistfile = open(outputprefix+"cdsexonendlist.txt","w")
    cdsexonendlistfile.write(cdsexonendlistfile_buffer)
    cdsexonendlistfile.close()

    segmentfile = open(outputprefix+"segment.txt","w")
    segmentfile.write(segmentfile_buffer )
    segmentfile.close()

    orthologygrouplistfile_buffer = ""
    if( extendedorthologygroups != []):
        for i in range(len(extendedorthologygroups)):
            orthologygrouplistfile_buffer += ">"+str(i)+"\n"
            for k in extendedorthologygroups[i]:
                cds = extendedsourcedata[k]
                cdsid,null,null ,null = cds
                orthologygrouplistfile_buffer += cdsid+"\n"
            orthologygrouplistfile_buffer += "\n"
            
        orthologygrouplistfile = open(outputprefix+"orthologygrouplist.txt","w")
        orthologygrouplistfile.write(orthologygrouplistfile_buffer)
        orthologygrouplistfile.close()

def pool_print_blocklist(data, outputformat, outputalignment):
    cdsid,geneid,cdsseq, cdsgeneid,geneseq, blocklist, status, cdsexon = data
    string_result_buffer = ""
    string_cdsexonendlist_buffer = ""
    segment_matches = []
    
    cdslength = len(cdsseq)
    string_result_buffer += cdsid + "\t" + geneid + "\t" + str(cdslength)  + "\t" + str("%.2f" % cover_percentage(blocklist,cdslength)) + "\t" + str(status)  + "\n"
    string_res, segments = print_blocklist(cdsid,cdsgeneid,geneid,cdsseq, geneseq, blocklist, outputformat, outputalignment,cdsexon)
    string_result_buffer += string_res
    segment_matches += segments
    
    if(cdsgeneid == geneid):
        string_cdsexonendlist_buffer += print_exonextremitylist(cdsid, geneid, blocklist)
    return [cdsid,geneid,string_result_buffer,string_cdsexonendlist_buffer, segment_matches]



def print_blocklist(cdsid, cdsgeneid,geneid, cds, gene, blocklist,outputformat,outputalignment,cdsexon):
    """
    This function is the main function of the graphic 
    representation of the prédiction.
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    blocklist:
    outfile: 
    outputformat:
    cdsexon:
        

    Returns
    -------
   
    
    """

    string_buffer = ""
    segment_matches = []
    cds_len = len(cds)

    if(len(blocklist) > 0):
        first_block = blocklist[0]
        first_qs = first_block[0]
        if(0 < first_qs):
            string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [0,first_qs], outputformat)
            string_buffer += string_to_print

        
        for i in range(len(blocklist)-1):
            block = blocklist[i]
            string_to_print,segments = compute_aln_string(cdsid, cdsgeneid, geneid, cds, gene, block, outputformat,outputalignment,cdsexon)
            string_buffer += string_to_print
            segment_matches += segments

            next_block = blocklist[i+1]
            block_qe = block[1]
            next_block_qs = next_block[0]
            if(block_qe < next_block_qs):
                string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [block_qe ,next_block_qs], outputformat)
                string_buffer += string_to_print

        last_block = blocklist[-1]
        string_to_print,segments = compute_aln_string(cdsid, cdsgeneid, geneid, cds, gene, last_block, outputformat,outputalignment,cdsexon)
        string_buffer += string_to_print
        segment_matches += segments

        last_block_qe = last_block[1]
        if(last_block_qe < cds_len):
            string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [last_block_qe,cds_len], outputformat)
            string_buffer += string_to_print
        
    else:
        string_to_print = compute_unaln_string(cdsid, geneid, cds, gene, [0 ,cds_len], outputformat)
        string_buffer += string_to_print
        
    #outfile.write("\n")
    return string_buffer, segment_matches

def print_wholealignment(cdsid, geneid, cds, gene, blocklist,outfile, outputformat):
    """
    This function
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    blocklist:
    outfile: 
    outputformat:
        

    Returns
    -------
   
    
    """
    cds_len = len(cds)
    gene_len = len(cds)
    cds_=""
    gene_ = ""

    sequence1 = ""
    sequence2 = ""
    if(len(cds_)==len(gene_)):
        sequence1 = gene_
        sequence2 = cds_
    elif(len(cds_)== 0):
        sequence1 = gene_
        sequence2 = '-' * len(sequence1)
    elif(len(gene_)== 0):
        sequence2 = cds_
        sequence1 = '-' * len(sequence2)
    else:
        alignment = pairwise2.align.globalms(gene_, cds_,2,0,-10,-1)
        sequence1, sequence2 = alignment[0][0],alignment[0][1]

        
    for i in range(len(blocklist)):
        block = blocklist[i]

        if (i==0):# first block
            cds_ += '-' * block[2]
            gene_ += gene[:block[2]]
            cds_ += cds[:block[0]]            
            gene_ += '-' * block[0]

        else:
            cds_ += '-' * (block[2] - blocklist[i-1][3])
            gene_ += gene[blocklist[i-1][3]:block[2]]
            cds_ += cds[blocklist[i-1][1]:block[0]]            
            gene_ += '-' * (block[0] - blocklist[i-1][1])
            
        sequence1 = ""
        sequence2 = ""
        cdsblock = cds[block[0]:block[1]]
        geneblock = gene[block[2]:block[3]]
        if(len(cdsblock)==len(geneblock)):
            sequence1 = geneblock
            sequence2 = cdsblock
        elif(len(cdsblock)== 0):
            sequence1 = geneblock
            sequence2 = '-' * len(sequence1)
        elif (len(geneblock)==0) :
            sequence2 = cdsblock
            sequence1 = '-' * len(sequence2)
        else:
            
            alignment = pairwise2.align.globalms(geneblock, cdsblock,2,0,-10,-1)
            sequence1, sequence2 = alignment[0][0],alignment[0][1]
        gene_ += sequence1
        cds_ += sequence2

        cds_ += '-' * (gene_len - blocklist[-1][3])
        gene_ += gene[blocklist[-1][3]:gene_len]
        cds_ += cds[blocklist[-1][1]:cds_len]            
        gene_ += '-' * (cds_len - blocklist[-1][1])

    outfile.write(">"+geneid+"\n")
    outfile.write(gene_+"\n")
    outfile.write(">"+cdsid+"\n")
    outfile.write(cds_+"\n")

          
def compute_aln_string(cdsid, cdsgeneid,geneid, cds, gene,block, outputformat,outputalignment,cdsexon):
    """
    This function produce the visual representation for each aligned block using a global alignment
     
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    block:
    cdsexon:
     
    outputformat:
        

    Returns
    -------
    string_to_print
    
    """
    string_to_print = ""
    
    block_qs = block[0] #query start
    block_qe = block[1] #query start
    block_ss = block[2] #subject start
    block_se = block[3] #subject end
    #block_identity = "%.2f" % (compute_block_identity(cds, gene,block))
    gene_= gene[block_ss:block_se]
    cds_= cds[block_qs:block_qe]

    sequence1 = ""
    sequence2 = ""
    block_identity = 0.0
    if(len(cds_)==len(gene_)):
        sequence1 = gene_
        sequence2 = cds_
    elif(len(cds_)== 0):
        sequence1 = gene_
        sequence2 = '-' * len(sequence1)
    elif(len(gene_)== 0):
        sequence2 = cds_
        sequence1 = '-' * len(sequence2)
    else:
        if(outputalignment == "zs"):
            alignment = pairwise2.align.globalms(gene_, cds_,2,0,-10,-1)
            sequence1, sequence2 = alignment[0][0],alignment[0][1]
        elif(outputalignment == "fsepsa"):
            #alignment = pairwise2.align.globalms(gene_, cds_,2,0,-10,-1)
            #sequence1, sequence2 = alignment[0][0],alignment[0][1]
            fsopen= -30
	    gapopen= -11
	    gapextend=-1
	    fsextend=-1
	    saa = ScoringMatrix('src/fsepsa/ressources/BLOSUM62.txt')
	    saa.load()
	    san = ScoringMatrix()
	    san.init_similarity()
	    arg = [fsopen, gapopen, gapextend, fsextend ]
	    score, sequence1, sequence2 = fse(gene_, cds_, arg, saa, san)

    aln_length = len(sequence1)

    block_identity = "%.2f" % (1.0 * computeAlignmentPercentIdentity(sequence1, sequence2) /100)

    if(cdsgeneid==geneid):
        assert(block_identity == "1.00")

    segment_matches = compute_segment_matches(sequence1, sequence2, block_ss, block_qs, block_identity)

    segment_matches = cut_segment_matches(segment_matches,cdsexon)
    
    if(outputformat == "list" or outputformat == "aln"):
        string_to_print = cdsid + "\t" + geneid + "\t" + str(aln_length) + "\t" + str(block_qs) + "\t" + str(block_qe) + "\t" + str(block_ss) +  "\t" + str(block_se) + "\t" + str(block_identity) +  "\t" + gene[block_ss-2:block_ss] + "<Exon>" + gene[block_se:block_se+2] + "\n"
    
    if(outputformat == "aln"):
        sequence1 = gene[block_ss-BORDER_LENGTH:block_ss] + sequence1 + gene[block_se:block_se+BORDER_LENGTH]
        sequence2 = BORDER_LENGTH*" " + sequence2 + BORDER_LENGTH*" "

        aln_srspair = format_alignment(sequence1,sequence2)

        string_to_print +=  aln_srspair
    if(outputformat == "listaln"):
        #string_to_print+=string(segment_matches)+"\n"
        #print segment_matches
        for match in segment_matches:
            genestart,cdsstart,length,pid = match
            string_to_print +=cdsid + "\t" + geneid + "\t" + str(length) + "\t" + str(cdsstart) + "\t" + str(cdsstart+length) + "\t" + str(genestart) +  "\t" + str(genestart+length) + "\t" + str(pid) +  "\n"
                    
            #string_to_print +=cdsid + "\t" + geneid + "\t" + str(length) + "\t" + str(cdsstart) + "\t" + str(cdsstart+length) + "\t" + str(genestart) +  "\t" + str(genestart+length) + "\t" + str(pid) +  "\n"
        
    return string_to_print, segment_matches

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
                
  
def compute_unaln_string(cdsid, geneid, cds, gene,interval, outputformat):
    """
    This function produce the visual representation for each unaligned cds block
       
    Parameters
    ----------

    
    cdsid:
    geneid:
    cds: 
    gene:
    interval:
     
    outputformat:
        

    Returns
    -------
    string_to_print
    
   
    """
    string_to_print = ""
    
    block_qs = interval[0] #query start
    block_qe = interval[1] #query end

    length = block_qe - block_qs
    
    string_to_print = cdsid + "\t" + geneid + "\t" + str(length) + "\t" + str(block_qs) + "\t" + str(block_qe) + "\t" + "-" +  "\t" + "-" + "\n"
    
    if(outputformat == "aln"):

        cds_= cds[block_qs:block_qe]
        sequence = format_sequence(cds_)

        string_to_print +=  sequence
        
    return string_to_print

def format_alignment(sequence1,sequence2):
    """
    This function is use to produce the visualization format with terminal adapt line length
       
    Parameters
    ----------
    sequence1:
    sequence2:
     
           

    Returns
    -------
    aln_srspair:
    
    
    """

    aln_srspair = ""

    markup_line = compute_markup_line(sequence1,sequence2)
    
    start = 0
    while start < len(sequence1):
        line_length = LENGTH_ALN_LINE
        if(start+line_length > len(sequence1)):
            line_length = len(sequence1)-start
        subsequence1 = sequence1[start:start+line_length]
        subsequence2 = sequence2[start:start+line_length]
        submarkup_line = markup_line[start:start+line_length]
        
        aln_srspair +=  subsequence1 + "\n"
        aln_srspair +=  submarkup_line + "\n"
        aln_srspair += subsequence2  + "\n\n"
        start += line_length

    return aln_srspair

def format_sequence(seq):
    """
    This function 

    Parameters
    ----------
    seq:
     
           

    Returns
    -------
    sequence:
    
    
    """
    sequence = ""
    
    start = 0
    while start < len(seq):
        line_length = LENGTH_ALN_LINE
        if(start+line_length > len(seq)):
            line_length = len(seq)-start
        subseq = seq[start:start+line_length]
        
        sequence +=  subseq + "\n"
        start += line_length

    sequence += "\n"
    return sequence


def compute_markup_line(sequence1,sequence2):
    """
    This function is use to return the central line of the visual alignment with the match and mismatch
    Parameters
    ----------
    sequence1:
    sequence2:
         

    Returns
    -------
        
    markup_line:
    
    """

    markup_line = ""
    for i in range(len(sequence1)):
        if(sequence1[i]==sequence2[i]):
            markup_line += "|"
        else:
            markup_line += " "
    return markup_line

def print_exonextremitylist(cdsid, geneid, blocklist):#, cdsexonendlistfile,geneexonstartlistfile,geneexonendlistfile):
    """
    This function 

    Parameters
    ----------
    cdsid:
    geneid:
    blocklist:
    cdsexonendlistfile:
    geneexonstartlistfile:
    geneexonendlistfile:
         

    Returns
    -------
        
    
    """

    string_buffer = ""
    string_buffer += ">"+cdsid+"\n"
    #geneexonstartlistfile.write(">"+geneid+"\n")
    #geneexonendlistfile.write(">"+geneid+"\n")
    for block in  blocklist:
        string_buffer += str(block[1])+" "
        #geneexonstartlistfile.write(str(block[2])+" ")
        #geneexonendlistfile.write(str(block[3])+" ")
    string_buffer += "\n"
    #geneexonstartlistfile.write("\n")
    #geneexonendlistfile.write("\n")
    return string_buffer
