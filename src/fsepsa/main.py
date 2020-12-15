#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os

from scoring_matrix import *

from fse import fse

from translator import aamap
def main():
	'''
	src/fse_main.py [-h] [-go GAPOPEN] [-ge GAPEXTEND] [-fso FSOPEN]
				          [-fse FSEXTEND] [-aa AMINOACIDMATRIX]
				          [-d DATASEQUENCE DATAFORMAT] [-o OUTFILE]  [-of OUTFORMAT]
	−30 for fs_open_
	cost and −1 for fs_extend_cost.
				         
	'''
	seq1='ACTGACACTAAGATAGACACAGATAGACACAGATAGAGAGAGACGTTTTGGGGTAAAGTCCCTGGTGTA'
	seq2='ACACACTAAGATAGACACAGATAGACACAGATAGAGAGAGACGTTTTGGGGTAAAGTCCCTGGTGTAAA'
	fsopen= -30
	gapopen= -11
	gapextend=-1
	fsextend=-1
	saa = ScoringMatrix('./ressources/BLOSUM62.txt')

	saa.load()

	san = ScoringMatrix()

	san.init_similarity()

	arg = [fsopen, gapopen, gapextend, fsextend ] 
	score, res_1, res_2 = fse(seq1, seq2, arg, saa, san)
	print('score, res_1, res_2', score, res_1, res_2)
main()
