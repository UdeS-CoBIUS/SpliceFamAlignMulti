#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``fse.py`` **module description**:

This module implements Theorem1 describing the dynamic programming algorithm in manuscript:
    * S. Jammali et al. Aligning protein-coding sequences with frameshift 
    extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

import numpy as np

from Lemma1 import lemma1
from Lemma2 import lemma2
from LemmaG import lemma_ga
from LemmaG import lemma_gb
from Case_D import *
from Case_Df import *
from Case_G import *

def gap(length, gap_open_cost, gap_extension_cost):
    value = 0
    if(length > 0):
        value = gap_open_cost + length * gap_extension_cost
    return value


def init_fse(seq_a, seq_b, fsopen, gap_open_cost, gap_extension_cost, sub_an):
    n = len(seq_a)
    m = len(seq_b)

    table_d[0, 0] = 0
    table_ga[0, 0] = -np.inf
    table_gb[0, 0] = -np.inf

    for i in range(1, n):
        table_d[i, 0] = gap(np.floor(i/3), gap_open_cost, gap_extension_cost)
        table_df[i, 0] = table_d[i, 0]

        if i % 3 == 1:
            table_df[i, 0] += sub_an[seq_a[i+1], seq_b[1]]/2 + sub_an[seq_a[i+2], seq_b[2]]/2 + fsopen
        elif i % 3 == 2:
            table_df[i, 0] += sub_an[seq_a[i+1], seq_b[1]]/2 + fsopen

        table_gb[i, 0] = -np.inf
        
    for j in range(1, m):
        table_d[0, j] = gap(np.floor(j/3), gap_open_cost, gap_extension_cost)
        table_df[0, j] = table_d[0, j]

        if j % 3 == 1:
            table_df[0, j] += sub_an[seq_a[1], seq_b[1+j]]/2 + sub_an[seq_a[2], seq_b[j+2]]/2 + fsopen
        elif j % 3 == 2:
            table_df[0, j] += sub_an[seq_a[1], seq_b[1+j]]/2 + fsopen

        table_ga[0, j] = -np.inf

def init_fct_vector():
    global d_fct
    global df_fct
    global ga_fct
    global gb_fct

    d_fct = [case_d1_1, case_d1_2, case_d1_3, case_d1_4, case_d1_5, case_d1_6, case_d1_7, case_d1_8, case_d1_9,
             case_d1_10, case_d1_11, case_d1_12, case_d1_13, case_d1_14, case_d1_15, case_d1_16, case_d1_17, case_d1_18,
             case_d2_1, case_d2_2, case_d2_3, case_d2_4, case_d2_5, case_d2_6, case_d2_7, case_d2_8, case_d2_9,
             case_d2_10, case_d2_11, case_d3_1, case_d3_2, case_d3_3, case_d3_4, case_d3_5, case_d3_6, case_d3_7,
             case_d3_8, case_d3_9, case_d3_10, case_d3_11, case_d4_1, case_d4_2, case_d4_3]

    df_fct = [case_df1_1, case_df2_1, case_df2_2, case_df2_3, case_df2_4, case_df2_5, case_df3_1, case_df3_2,
              case_df3_3, case_df3_4, case_df3_5, case_df4_1, case_df4_2, case_df4_3, case_df5_1, case_df5_2,
              case_df5_3]

    ga_fct = [case_ga_1, case_ga_2]

    gb_fct = [case_gb_1, case_gb_2]

def backtrack(seq_a, seq_b):
    i = len(seq_a) - 1
    j = len(seq_b) - 1
    res_1 = res_2 = ''
    memories = [memory, memory_df, memory_ga, memory_gb]
    fct_vect = [d_fct, df_fct, ga_fct, gb_fct]
    ind = 0
    init_fct_vector()

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            idx = memories[ind][i, j]
            res_1, res_2, i, j, ind = fct_vect[ind][idx](seq_a, seq_b, res_1, res_2, i, j)

        elif i > 0:
            res_1 = seq_a[i] + res_1
            res_2 = '-' + res_2
            i -= 1

        elif j > 0:
            res_1 = '-' + res_1
            res_2 = seq_b[j] + res_2
            j -= 1

    return res_1, res_2


def fse(seq_a, seq_b, arg, sub_aa, sub_an):
    [fsopen, gapopen, gapextend, fsextend ]= arg 

    global table_d
    global table_df
    global table_ga
    global table_gb
    global memory
    global memory_df
    global memory_ga
    global memory_gb

    n = len(seq_a)
    m = len(seq_b)
    seq_a = ' ' + seq_a
    seq_b = ' ' + seq_b
    table_d = np.zeros((n+1, m+1))
    table_df = np.zeros(table_d.shape)
    table_ga = np.zeros(table_d.shape)
    table_gb = np.zeros(table_d.shape)
    memory = np.zeros(table_d.shape, dtype='int')
    memory_df = np.zeros(table_d.shape, dtype='int')
    memory_ga = np.zeros(table_d.shape, dtype='int')
    memory_gb = np.zeros(table_d.shape, dtype='int')
    
    # for j in xrange(1, m+1):
    #     table_ga[0][j] = -np.inf
    # for i in xrange(1, n+1):
    #     table_gb[i][0] = -np.inf
        
    init_fse(seq_a, seq_b, fsopen, gapopen, gapextend, sub_an)

    for i in range(1, n+1):
        for j in range(1, m+1):
            lemma_ga(seq_a, seq_b, gapopen, gapextend, table_ga, table_d, memory_ga, i, j)
            lemma_gb(seq_a, seq_b, gapopen, gapextend, table_gb, table_d, memory_gb, i, j)
            lemma1(seq_a, seq_b, fsopen, fsextend, sub_aa, sub_an, table_d,
                   table_df, table_ga, table_gb, memory, i, j)
            lemma2(seq_a, seq_b, fsopen, fsextend, sub_aa, sub_an, table_d, table_df, memory_df, i, j)

    init_fct_vector()
    score = table_d[n, m]
    res_1, res_2 = backtrack(seq_a, seq_b)    

    return score, res_1, res_2
