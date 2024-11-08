#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``translator.py`` **module description**:

This module implements the translation of coding sequences into amino acid sequences.

From http://thepuri-st.blogspot.ca/2013/03/bioinformatics-adventures-in-python.html

.. moduleauthor:: DrG

"""

# Create a dictionary of amino acid codes
aamap = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}


def prot(str):
    # Decode RNA into protein sequence

    # FIXME retourner une etoile a la place de retourner un 'stop'

    # Read in file
    f = open(handle, 'r')
    mRNA = f.readline()

    protein = []
    for codon in rnareader(mRNA):

        # if we hit a stop codon
        if aamap[codon] == "*":
            break

        if codon == "AUG":
            start = True

        if start:
            # if we hit a stop codon
            if aamap[codon] == "*":
                break
            aa = aamap[codon]
            protein.append(aa)


    print(''.join(protein))

# Also need to build a generator to return
def rnareader(rnastr):
    aacodon = []
    for nuc in rnastr:
        aacodon.append(nuc)
        if len(aacodon) >= 3:
            yield ''.join(aacodon)
            aacodon = []
