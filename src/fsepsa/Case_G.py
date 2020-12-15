#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``Case_G.py`` **module description**:

This module implements the backtracking steps for dynamic programming tables GA and GB Lemma2 used to handle affine gap costs in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: Aïda Ouangraoua

February 2016

"""

def case_ga_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 2

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '---' + res_b

    i -= 3

    return res_a, res_b, i, j, ind

def case_ga_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '---' + res_b

    i -= 3

    return res_a, res_b, i, j, ind

def case_gb_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 3

    res_a = '---' + res_a
    res_b = seq_b[j-2:j+1] + res_b

    j -= 3

    return res_a, res_b, i, j, ind

def case_gb_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '---' + res_a
    res_b = seq_b[j-2:j+1] + res_b

    j -= 3

    return res_a, res_b, i, j, ind

