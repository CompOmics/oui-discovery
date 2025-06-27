#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
List of enzymes and the regex that can be selected to perform the in silico
digestion.
"""
# ------------------------------------------------------------------------------

# Regular expression used to find the cleavage site
cleavage_site = {
            'trypsin-p': '(?<=[RK](?!P))',      # Cleaves after arginine or lysine not followed by proline
            'trypsin': '(?<=[RK])',      # Cleaves after arginine or lysine
            'lysc-p' : '(?<=K(?!P))',   # Cleaves after lysine not followed by proline
            'lysc': '(?<=K)',    # Cleaves after lysine
            'lysn': '(?=K)',      # Cleaves before lysine
            'argc': '(?<=R)',    # Cleaves after arginine
            'aspc': '(?<=D)',   # Cleaves after aspartic acid
            'aspn': '(?=D)',      # Cleaves before aspartic acid
            'gluc': '(?<=E)',    # Cleaves after glutamic acid
            'glun': '(?=E)',      # Cleaves before glutamic acid
            'chymotrypsin+': '(?<=[YWFLM])',     # Cleaves after tyrosine, tryptophane, phenylalanine, leucine, methionine
            'chymotrypsin': '(?<=[YWF])'     # Cleaves after tyrosine, tryptophane, phenylalanine
        }
