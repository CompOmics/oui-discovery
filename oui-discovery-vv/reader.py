#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
Some functions to be used for the digestion 
"""
# ------------------------------------------------------------------------------
from collections import defaultdict
# ------------------------------------------------------------------------------

def fasta_to_dict(fasta_file):
    """ Convert fasta file to python dictionary """
    fasta_dict = defaultdict(str)
    with open(fasta_file, 'r', encoding='UTF-8') as file:
        for line in file:
            if line.startswith('>'):
                id_ = line.split(' ')[0].strip('>')
            else:
                fasta_dict[id_] += line.strip('\n')
    return fasta_dict
