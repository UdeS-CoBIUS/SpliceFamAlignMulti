#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``Case_Df.py`` **module description**:

This module implements the backtracking steps for Lemma2 used to fill dynamic programming table Df of the algorithm described in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

def case_df1_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    return res_a, res_b, i, j, ind


def case_df2_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_a = seq_a[i-1:i+1] + res_a
    res_b = seq_b[j-1:j+1] + res_b
    
    i -= 2
    j -= 2
    
    return res_a, res_b, i, j, ind


def case_df2_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_a = seq_a[i-1:i+1] + res_a
    res_b = '-' + seq_b[j] + res_b
    
    i -= 2 
    j -= 1
    
    return res_a, res_b, i, j, ind


def case_df2_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_a = seq_a[i-1:i+1] + res_a
    res_b = seq_b[j] + '-' + res_b
    
    i -= 2
    j -= 1
    
    return res_a, res_b, i, j, ind


def case_df2_4(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_a = seq_a[i-1:i+1] + res_a
    res_b = '--' + res_b
    i -= 2
    
    return res_a, res_b, i, j, ind


def case_df2_5(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    return res_a, res_b, i, j, ind


def case_df3_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_b = seq_b[j-1:j+1] + res_b
    res_a = seq_a[i-1:i+1] + res_a
    
    j -= 2
    i -= 2
    
    return res_a, res_b, i, j, ind


def case_df3_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_b = seq_b[j-1:j+1] + res_b
    res_a = '-' + seq_a[i] + res_a
    
    j -= 2 
    i -= 1
    
    return res_a, res_b, i, j, ind


def case_df3_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_b = seq_b[j-1:j+1] + res_b
    res_a = seq_a[i] + '-' + res_a
    
    j -= 2
    i -= 1
    
    return res_a, res_b, i, j, ind


def case_df3_4(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_b = seq_b[j-1:j+1] + res_b
    res_a = '--' + res_a
    j -= 2
    
    return res_a, res_b, i, j, ind


def case_df3_5(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    return res_a, res_b, i, j, ind


def case_df4_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_a = seq_a[i] + res_a
    res_b = seq_b[j] + res_b
    
    i -= 1
    j -= 1
    return res_a, res_b, i, j, ind


def case_df4_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_a = seq_a[i] + res_a
    res_b = '-' + res_b
    
    i -= 1
    
    return res_a, res_b, i, j, ind


def case_df4_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    return res_a, res_b, i, j, ind


def case_df5_1(seq_a, seq_b, res_a, res_b, i, j):
    ind = 1
    
    res_b = seq_b[j] + res_b
    res_a = seq_a[i] + res_a
    
    j -= 1
    i -= 1
    return res_a, res_b, i, j, ind


def case_df5_2(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    
    res_b = seq_b[j] + res_b
    res_a = '-' + res_a
    
    j -= 1
    
    return res_a, res_b, i, j, ind


def case_df5_3(seq_a, seq_b, res_a, res_b, i, j):
    ind = 0
    return res_a, res_b, i, j, ind
