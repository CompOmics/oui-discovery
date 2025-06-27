#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
Dict of the molecular mass (in dalton) of each amino acid.
Reference: matrixscience.com/help/aa_help.html
"""
# ------------------------------------------------------------------------------
atom = {
            'H': 1.00782503521,
            'O': 15.9949146221,
            'C': 12.0000000,
            'N': 14.0030740052,
            'P': 30.97376151,
            'S': 31.97207069,
            'Se': 79.916522,
            }

aminoAcid = {
            "A": atom['C']*3 + atom['H']*5 + atom['N'] + atom['O'],
            "R": atom['C']*6 + atom['H']*12 + atom['N']*4 + atom['O'],
            "N": atom['C']*4 + atom['H']*6 + atom['N']*2 + atom['O']*2,
            "D": atom['C']*4 + atom['H']*5 + atom['N'] + atom['O']*3,
            "C": atom['C']*3 + atom['H']*5 + atom['N'] + atom['O'] + atom['S'],
            "E": atom['C']*5 + atom['H']*7 + atom['N'] + atom['O']*3,
            "Q": atom['C'] *5+ atom['H']*8 + atom['N']*2 + atom['O']*2,
            "G": atom['C']*2 + atom['H']*3 + atom['N'] + atom['O'],
            "H": atom['C']*6 + atom['H']*7 + atom['N']*3 + atom['O'],
            "I": atom['C']*6 + atom['H']*11 + atom['N'] + atom['O'],
            "L": atom['C']*6 + atom['H']*11 + atom['N'] + atom['O'],
            "K": atom['C']*6 + atom['H']*12 + atom['N']*2 + atom['O'],
            "M": atom['C']*5 + atom['H']*9 + atom['N'] + atom['O'] + atom['S'],
            "F": atom['C']*9 + atom['H']*9 + atom['N'] + atom['O'],
            "P": atom['C']*5 + atom['H']*7 + atom['N'] + atom['O'],
            "S": atom['C']*3 + atom['H']*5 + atom['N'] + atom['O']*2,
            "T": atom['C']*4 + atom['H']*7 + atom['N'] + atom['O']*2,
            'U': atom['C']*3 + atom['H']*7 + atom['N'] + atom['O']*2 + atom['Se'],
            "W": atom['C']*11 + atom['H']*10 + atom['N']*2 + atom['O'],
            "Y": atom['C']*9 + atom['H']*9 + atom['N'] + atom['O']*2,
            'X': 120.1779,  # Mean Mass #@uncomment
            "V": atom['C']*5 + atom['H']*9 + atom['N'] + atom['O'],
            "Z": atom['C']*5 + atom['H']*7 + atom['N'] + atom['O']*3, #@ glutamic acid or glutamine, so put Z=E 
            "B": atom['C']*4 + atom['H']*6 + atom['N']*2 + atom['O']*2, #@Aspartic acid or Asparagine B=N
            "J": atom['C']*6 + atom['H']*11 + atom['N'] + atom['O'] #@Leucine or Isoleucine J=L
        }
