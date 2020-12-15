#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``Lemma1.py`` **module description**:

This module implements Lemma1 used to fill dynamic programming table D 
of the algorithm described in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

.. moduleauthor:: François Bélanger

February 2016

"""

import numpy as np
from translator import aamap

def lemma1(seq_a, seq_b, fs_open_cost, fs_extension_cost, sub_aa, sub_an, table_d, table_df, table_ga, table_gb, memory, i, j):
    if i % 3 == 0 and j % 3 == 0:
        values = np.zeros(18)
        values[:] = -np.inf

        values[0] = sub_aa[aamap[str(seq_a[i-2:i+1])], aamap[str(seq_b[j-2:j+1])]] + table_d[i-3, j-3]
        values[1] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-1]] + table_d[i-3, j-2] + 2 * fs_open_cost
        values[2] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-2], seq_b[j-1]] + table_d[i-3, j-2] + 2 * fs_open_cost
        values[3] = sub_an[seq_a[i], seq_b[j]] + table_d[i-3, j-1] + 2 * fs_open_cost
        values[4] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-1]] + table_d[i-2, j-3] + 2 * fs_open_cost
        values[5] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-2]] + table_d[i-2, j-3] + 2 * fs_open_cost
        values[6] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-3] + 2 * fs_open_cost
        values[7] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-1] + 2 * fs_open_cost

        values[8] = sub_an[seq_a[i-1], seq_b[j]]/2 + sub_an[seq_a[i-2], seq_b[j-1]]/2 + table_df[i-3, j-2] + \
            fs_open_cost

        values[9] = sub_an[seq_a[i-1], seq_b[j]] + table_d[i-3, j-1] + 2 * fs_open_cost
        values[10] = sub_an[seq_a[i-2], seq_b[j]]/2 + table_df[i-3, j-1] + fs_open_cost
        values[11] = table_ga[i, j]
        values[12] = table_d[i-1, j] + fs_open_cost

        values[13] = sub_an[seq_a[i], seq_b[j-1]]/2 + sub_an[seq_a[i-1], seq_b[j-2]]/2 + table_df[i-2, j-3] + \
            fs_open_cost

        values[14] = sub_an[seq_a[i], seq_b[j-1]] + table_d[i-1, j-3] + 2 * fs_open_cost
        values[15] = sub_an[seq_a[i], seq_b[j-2]]/2 + table_df[i-1, j-3] + fs_open_cost
        values[16] = table_gb[i, j]
        values[17] = table_d[i, j-1] + fs_open_cost

        memory[i, j] = np.argmax(values)
        table_d[i, j] = np.amax(values)


    elif i % 3 == 0 and j % 3 != 0:
        values = np.zeros(11)
        values[:] = -np.inf

        if i >= 3 and j >= 3:
            values[0] = sub_aa[aamap[str(seq_a[i-2:i+1])], aamap[str(seq_b[j-2:j+1])]]/2 + table_df[i-3, j-3] + \
                fs_extension_cost + sub_an[seq_a[i], seq_b[j]]/2

            if (j - 1) % 3 != 0:
                values[0] += sub_an[seq_a[i-1], seq_b[j-1]]/2

        if i >= 3 and j >= 2:
            open_mult = 2 if (j - 1) % 3 == 0 else 1
            values[1] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-1]] + table_d[i-3, j-2] + \
                open_mult * fs_open_cost

            values[2] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-2], seq_b[j-1]] + table_df[i-3, j-2] + fs_open_cost
            if (j - 1) % 3 == 0:
                values[2] -= (sub_an[seq_a[i-2], seq_b[j-1]]/2)

        if j >= 1:
            if i >= 3:
                values[3] = sub_an[seq_a[i], seq_b[j]] + table_d[i-3, j-1] + fs_open_cost

            values[4] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-1] + fs_open_cost

        if i >= 3 and j >= 2:
            values[5] = sub_an[seq_a[i-1], seq_b[j]] + sub_an[seq_a[i-2], seq_b[j-1]] + table_df[i-3, j-2] + \
                fs_open_cost

            if (j - 1) % 3 == 0:
                values[5] -= (sub_an[seq_a[i-2], seq_b[j-1]]/2)

        if i >= 3:
            values[6] = sub_an[seq_a[i-1], seq_b[j]] + table_d[i-3, j-1] + fs_open_cost
            values[7] = sub_an[seq_a[i-2], seq_b[j]] + table_d[i-3, j-1] + fs_open_cost

        if i >= 3:
            values[8] = table_ga[i, j]

        values[9] = table_d[i-1, j] + fs_open_cost
        values[10] = table_d[i, j-1]

        memory[i, j] = np.argmax(values) + 18
        table_d[i, j] = np.max(values)

    elif i % 3 != 0 and j % 3 == 0:
        values = np.zeros(11)
        values[:] = -np.inf

        if i >= 3 and j >= 3:
            values[0] = sub_aa[aamap[str(seq_a[i-2:i+1])], aamap[str(seq_b[j-2:j+1])]]/2 + table_df[i-3, j-3] + \
                fs_extension_cost + sub_an[seq_a[i], seq_b[j]]/2

            if (i - 1) % 3 != 0:
                values[0] += sub_an[seq_a[i-1], seq_b[j-1]]/2

        open_mult = 2 if (i - 1) % 3 == 0 else 1

        if i >= 2 and j >= 3:
            values[1] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-1]] + table_d[i-2, j-3] + \
                open_mult * fs_open_cost

            values[2] = sub_an[seq_a[i], seq_b[j]] + sub_an[seq_a[i-1], seq_b[j-2]] + table_df[i-2, j-3] + fs_open_cost

            if (i - 1) % 3 == 0:
                values[2] -= (sub_an[seq_a[i-1], seq_b[j-2]]/2)

        if j >= 3:
            values[3] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-3] + fs_open_cost

        if j >= 1:
            values[4] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-1] + fs_open_cost

        if i >= 2 and j >= 3:
            values[5] = sub_an[seq_a[i], seq_b[j-1]] + sub_an[seq_a[i-1], seq_b[j-2]] + table_df[i-2, j-3] + \
                fs_open_cost

            if (i - 1) % 3 == 0:
                values[5] -= (sub_an[seq_a[i-1], seq_b[j-2]]/2)

        if j >= 3:
            values[6] = sub_an[seq_a[i], seq_b[j-1]] + table_d[i-1, j-3] + fs_open_cost
            values[7] = sub_an[seq_a[i], seq_b[j-2]] + table_d[i-1, j-3] + fs_open_cost

            values[8] = table_gb[i, j]

        values[9] = table_d[i, j-1] + fs_open_cost
        values[10] = table_d[i-1, j]

        memory[i, j] = np.argmax(values) + 18 + 11
        table_d[i, j] = np.max(values)

    elif i % 3 != 0 and j % 3 != 0:
        values = np.zeros(3)
        values[:] = -np.inf

        values[0] = sub_an[seq_a[i], seq_b[j]] + table_d[i-1, j-1]
        values[1] = table_d[i-1, j]
        values[2] = table_d[i, j-1]

        memory[i, j] = np.argmax(values) + 18 + 11 + 11
        table_d[i, j] = np.amax(values)

    return table_d
