#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``Case_D.py`` **module description**:

This module implements the backtracking steps for Lemma1 used to fill dynamic programming table D of the algorithm described in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

def case_d1_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 3
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '-' + seq_b[j-1:j+1] + res_b

    i -= 3
    j -= 2
    return res_a, res_b, i, j, ind


def case_d1_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-1] + '-' + seq_b[j] + res_b

    i -= 3
    j -= 2

    return res_a, res_b, i, j, ind


def case_d1_4(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '--' + seq_b[j] + res_b

    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d1_5(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '-' + seq_a[i-1:i+1] + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 2
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_6(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-1] + '-' + seq_a[i] + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 2
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_7(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '--' + seq_a[i] + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 1
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_8(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = seq_b[j] + res_b

    i -= 1
    j -= 1

    return res_a, res_b, i, j, ind


def case_d1_9(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-1:j+1] + '-' + res_b

    i -= 3
    j -= 2

    return res_a, res_b, i, j, ind


def case_d1_10(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '-' + seq_b[j] + '-' + res_b

    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d1_11(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j] + '--' + res_b

    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d1_12(seq_a, seq_b, res_a, res_b, i, j):
    ind = 2

    return res_a, res_b, i, j, ind


def case_d1_13(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = '-' + res_b

    i -= 1

    return res_a, res_b, i, j, ind


def case_d1_14(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-1:i+1] + '-' + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 2
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_15(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '-' + seq_a[i] + '-' + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 1
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_16(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i] + '--' + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 1
    j -= 3

    return res_a, res_b, i, j, ind


def case_d1_17(seq_a, seq_b, res_a, res_b, i, j):
    ind = 3

    return res_a, res_b, i, j, ind


def case_d1_18(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '-' + res_a
    res_b = seq_b[j] + res_b

    j -= 1

    return res_a, res_b, i, j, ind


def case_d2_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-2:j+1] + res_b

    i -= 3
    j -= 3

    return res_a, res_b, i, j, ind


def case_d2_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '-' + seq_b[j-1:j+1] + res_b
    
    i -= 3
    j -= 2

    return res_a, res_b, i, j, ind


def case_d2_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-1] + '-' + seq_b[j] + res_b
    
    i -= 3
    j -= 2

    return res_a, res_b, i, j, ind


def case_d2_4(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '--' + seq_b[j] + res_b
    
    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d2_5(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = seq_b[j] + res_b
    
    i -= 1
    j -= 1

    return res_a, res_b, i, j, ind


def case_d2_6(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j-1:j+1] + '-' + res_b
    
    i -= 3
    j -= 2

    return res_a, res_b, i, j, ind


def case_d2_7(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = '-' + seq_b[j] + '-' + res_b
    
    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d2_8(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i-2:i+1] + res_a
    res_b = seq_b[j] + '--' + res_b

    i -= 3
    j -= 1

    return res_a, res_b, i, j, ind


def case_d2_9(seq_a, seq_b, res_a, res_b, i, j):
    ind = 2

    return res_a, res_b, i, j, ind


def case_d2_10(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = '-' + res_b

    i -= 1

    return res_a, res_b, i, j, ind


def case_d2_11(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = '-' + res_a
    res_b = seq_b[j] + res_b

    j -= 1

    return res_a, res_b, i, j, ind


def case_d3_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_b = seq_b[j-2:j+1] + res_b
    res_a = seq_a[i-2:i+1] + res_a

    j -= 3
    i -= 3

    return res_a, res_b, i, j, ind


def case_d3_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j-2:j+1] + res_b
    res_a = '-' + seq_a[i-1:i+1] + res_a
    
    j -= 3
    i -= 2

    return res_a, res_b, i, j, ind


def case_d3_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_b = seq_b[j-2:j+1] + res_b
    res_a = seq_a[i-1] + '-' + seq_a[i] + res_a
    
    j -= 3
    i -= 2

    return res_a, res_b, i, j, ind


def case_d3_4(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j-2:j+1] + res_b
    res_a = '--' + seq_a[i] + res_a
    
    j -= 3
    i -= 1

    return res_a, res_b, i, j, ind


def case_d3_5(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j] + res_b
    res_a = seq_a[i] + res_a
    
    j -= 1
    i -= 1

    return res_a, res_b, i, j, ind


def case_d3_6(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1

    res_b = seq_b[j-2:j+1] + res_b
    res_a = seq_a[i-1:i+1] + '-' + res_a
    
    j -= 3
    i -= 2

    return res_a, res_b, i, j, ind


def case_d3_7(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j-2:j+1] + res_b
    res_a = '-' + seq_a[i] + '-' + res_a
    
    j -= 3
    i -= 1

    return res_a, res_b, i, j, ind


def case_d3_8(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j-2:j+1] + res_b
    res_a = seq_a[i] + '--' + res_a

    j -= 3
    i -= 1

    return res_a, res_b, i, j, ind


def case_d3_9(seq_a, seq_b, res_a, res_b, i, j):
    ind = 3

    return res_a, res_b, i, j, ind


def case_d3_10(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j] + res_b
    res_a = '-' + res_a

    j -= 1

    return res_a, res_b, i, j, ind


def case_d3_11(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = '-' + res_b
    res_a = seq_a[i] + res_a

    i -= 1

    return res_a, res_b, i, j, ind


def case_d4_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = seq_b[j] + res_b

    i -= 1
    j -= 1

    return res_a, res_b, i, j, ind


def case_d4_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_a = seq_a[i] + res_a
    res_b = '-' + res_b
    i -= 1

    return res_a, res_b, i, j, ind


def case_d4_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0

    res_b = seq_b[j] + res_b
    res_a = '-' + res_a
    j -= 1

    return res_a, res_b, i, j, ind
