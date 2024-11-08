#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``LemmaG.py`` **module description**:

This module implements affine gap costs usage in the algorithm described in manuscript:
    * F. Bélanger, S. Jammali, A. Rachati, A. Ouangraoua. Aligning coding sequences with frameshift extension penalties. (2016).

Affine gap costs are handled by adding two dynamic programming matrices GA and GB.

.. moduleauthor:: Aïda Ouangraoua

February 2016

"""

import numpy as np
from translator import aamap

def lemma_ga(seq_a, seq_b, gap_open_cost, gap_extension_cost, table_ga, table_d, memory_ga, i, j):
    values = np.zeros(2)
    values[:] = -np.inf
    
    values[0] = table_ga[i-3, j] + gap_extension_cost
    values[1] = table_d[i-3, j] + gap_open_cost + gap_extension_cost
    
    memory_ga[i, j] = np.argmax(values)
    table_ga[i, j] = np.amax(values)

    return table_ga

def lemma_gb(seq_a, seq_b, gap_open_cost, gap_extension_cost, table_gb, table_d, memory_gb, i, j):
    values = np.zeros(2)
    values[:] = -np.inf
    
    values[0] = table_gb[i, j-3] + gap_extension_cost
    values[1] = table_d[i, j-3] + gap_open_cost + gap_extension_cost
    
    memory_gb[i, j] = np.argmax(values)
    table_gb[i, j] = np.amax(values)

    return table_gb
