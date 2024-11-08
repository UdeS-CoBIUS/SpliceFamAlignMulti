#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``scoring_matrix.py`` **module description**:

This module implements a class for subtitution matrices used in the dynamic programming algorithm described in manuscript:
    * F. Bélanger, A. Rachati, A. Ouangraoua. Aligning protein-coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

import numpy as np

class ScoringMatrix:

    def __init__(self, filename=None):
        self._filename = filename
        self._alphabet = None
        self._scoring_matrix = None

    def _read_matrix(self, scoring_file):
        if self._scoring_matrix is None:
            sz = len(self._alphabet)
            self._scoring_matrix = np.zeros((sz, sz))
            
        i = 0

        for line in scoring_file:
            values = line.split()
            self._scoring_matrix[i, :] = list(map(float, values[1:]))
	    
            i += 1

    def _build_alphabet(self, line):
	
        self._alphabet = {}
	
        keys = line.split()

        i = 0
	
        for k in keys:
            self._alphabet[k] = i
	   
            i += 1

    @staticmethod
    def _skip_comment(scoring_file):
        line = scoring_file.readline()

        while line[0] == '#':
            line = scoring_file.readline()

        return line

    def load(self, filename=None):
        must_load = False

        if filename is not None:
            must_load = True
            self._filename = filename
            self._alphabet = None

        if self._scoring_matrix is None or must_load is True:
            scoring_file = open(self._filename, 'r')

            line = self._skip_comment(scoring_file)
            self._build_alphabet(line)
            self._read_matrix(scoring_file)
            scoring_file.close()

    def init_similarity(self):
        self._alphabet = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self._scoring_matrix = np.ones((4, 4))
        self._scoring_matrix *= -1.0

        for i in range(4):
            self._scoring_matrix[i, i] = 1.0

    def __getitem__(self, item):
        return self._scoring_matrix[self._alphabet[item[0]], self._alphabet[item[1]]]
