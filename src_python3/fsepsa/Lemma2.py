#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``Lemma2.py`` **module description**:

This module implements Lemma2 used to fill dynamic programming table Df 
of the algorithm described in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

import numpy as np
from translator import aamap


def lemma2(seq_a, seq_b, fs_open_cost, fs_extension_cost, sub_aa, sub_an, table_d, table_df, memory_df, i, j):
    len_a = len(seq_a)
    len_b = len(seq_b)

    if i % 3 == 0 and j % 3 == 0:
        memory_df[i, j] = 0
        table_df[i, j] = table_d[i, j]

    elif i % 3 == 2 and j % 3 == 0:
        values = np.zeros(5)
        values[:] = -np.inf

        if 2 <= i < len_a-1 and 2 <= j < len_b-1:
            values[0] = sub_aa[aamap[str(seq_a[i-1:i+2])], aamap[str(seq_b[j-1:j+2])]]/2 + table_df[i-2, j-2] + fs_extension_cost

        if 2 <= i < len_a-1 and j < len_b-1:
            values[1] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + sub_an[seq_a[i], seq_b[j]] + table_d[i-2, j-1] + \
                2 * fs_open_cost

            values[2] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + sub_an[seq_a[i-1], seq_b[j]]/2 + table_df[i-2, j-1] + \
                fs_open_cost

            values[3] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i-2, j] + fs_open_cost

        if i < len_a-1 and j < len_b-1:
            values[4] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j] + fs_open_cost

        memory_df[i, j] = np.argmax(values) + 1
        table_df[i, j] = np.amax(values)

    elif i % 3 == 0 and j % 3 == 2:
        values = np.zeros(5)
        values[:] = -np.inf

        if 2 <= i < len_a-1 and 2 <= j < len_b-1:
            values[0] = sub_aa[aamap[str(seq_a[i-1:i+2])], aamap[str(seq_b[j-1:j+2])]]/2 + table_df[i-2, j-2] + fs_extension_cost

        if 2 <= j < len_b-1 and i < len_a-1:
            values[1] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-2] + \
                2 * fs_open_cost

            values[2] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + sub_an[seq_a[i], seq_b[j-1]]/2 + table_df[i-1, j-2] + \
                fs_open_cost

            values[3] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j-2] + fs_open_cost

        if i < len_a-1 and j < len_b-1:
            values[4] = sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j] + fs_open_cost

        memory_df[i, j] = np.argmax(values) + 1 + 5
        table_df[i, j] = np.amax(values)

    elif i % 3 == 1 and j % 3 == 0:
        values = np.zeros(3)
        values[:] = -np.inf

        if i < len_a-2 and j < len_b-2:
            values[0] = sub_aa[aamap[str(seq_a[i:i+3])], aamap[str(seq_b[j:j+3])]]/2 + table_df[i-1, j-1] + fs_extension_cost

        if i < len_a-2 and j < len_b-2:
            values[1] = sub_an[seq_a[i+2], seq_b[j+2]]/2 + sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i-1, j] + \
                fs_open_cost
            values[2] = sub_an[seq_a[i+2], seq_b[j+2]]/2 + sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j] + \
                fs_open_cost

        memory_df[i, j] = np.argmax(values) + 1 + 5 + 5
        table_df[i, j] = np.amax(values)

    elif i % 3 == 0 and j % 3 == 1:
        values = np.zeros(3)
        values[:] = -np.inf

        if i < len_a-2 and j < len_b-2:
            values[0] = sub_aa[aamap[str(seq_a[i:i+3])], aamap[str(seq_b[j:j+3])]]/2 + table_df[i-1, j-1] + fs_extension_cost

        if i < len_a-2 and j < len_b-2:
            values[1] = sub_an[seq_a[i+2], seq_b[j+2]]/2 + sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j-1] + \
                fs_open_cost
            values[2] = sub_an[seq_a[i+2], seq_b[j+2]]/2 + sub_an[seq_a[i+1], seq_b[j+1]]/2 + table_d[i, j] + \
                fs_open_cost

        memory_df[i, j] = np.argmax(values) + 1 + 5 + 5 + 3
        table_df[i, j] = np.amax(values)

    return table_df
