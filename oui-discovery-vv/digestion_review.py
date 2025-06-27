#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
This script was adapted from https://github.com/Francis-B/Kiwi.

This script contains functions and a class to perform an in silico digestion
and check the uniqueness of peptide sequences. 
"""
# ------------------------------------------------------------------------------
import re
import numpy as np
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
            'X': 120.1779,  # Mean Mass
            "V": atom['C']*5 + atom['H']*9 + atom['N'] + atom['O'],
            "Z": atom['C']*5 + atom['H']*7 + atom['N'] + atom['O']*3, #@ glutamic acid or glutamine, so put Z=E 
            "B": atom['C']*4 + atom['H']*6 + atom['N']*2 + atom['O']*2, #@Aspartic acid or Asparagine B=N
            "J": atom['C']*6 + atom['H']*11 + atom['N'] + atom['O'] #@Leucine or Isoleucine J=L
        }


class Digestion:
    """ Instantiate this class to perform in silico digestion and find unique
        peptide """

    def __init__(self, protein_sequences, min_length=7, max_length=None, max_miscleavages=1, max_mass=4600):
        self.params = {
                       'min_length': min_length,
                       'max_length': max_length,
                       'max_miscleavages': max_miscleavages,
                       'max_mass': max_mass,
                       }

        self.proteins = protein_sequences
        self.peptides = {id: [] for id in self.proteins}
        self.outdir = 'digested_peptides.csv'
        self.is_unique = {}


    def get_peptides_list(self):
        """" Return flat list of all peptides """
        return np.unique([pep for peptides in self.peptides.values() for pep in peptides])


    def _get_mass(self, seq):
        """Sum the mass of all residues and add the mass of H (N-terminal) and OH 
         (C-terminus) to obtain the sequence mass"""
        _mass = sum((aminoAcid[aa] for aa in seq)) + atom['H']*2 + atom['O']
        return round(_mass, 4)


    def _join_sequences(self, base_sequences, miss):
        """ Iterate through base sequences (fully digested peptide) and join miss+1 
            adjacent sequences """
        num_sequences = len(base_sequences)-miss  # Number of sequences to loop on
        return [''.join(base_sequences[n:n+miss+1]) for n in range(num_sequences)]


    def _is_good_peptide(self, seq):
        """ Filter sequences with given threshold(s)"""
        if self.params['min_length'] is not None and \
           self.params['min_length'] > len(seq):
            return False
        if self.params['max_length'] is not None and \
           self.params['max_length'] < len(seq):
            return False
        if self.params['max_mass'] is not None and \
           self.params['max_mass'] < self._get_mass(seq):
            return False

        return True


    def cleave_proteins(self):
        """ Perform an in silico digestion of proteins with previously selected
        or default parameters """
        regexp = '(?<=[RK])' # Trypsin cleavage site

        for id_, protein_sequences in self.proteins.items():
            # split protein into perfectly digested sequences
            base_sequences = re.split(regexp, protein_sequences)
            miss = 0
            while miss <= self.params['max_miscleavages']:
                # Join base sequences in function of actual miscleavage value
                all_peptides = self._join_sequences(base_sequences, miss)
                # Filter all peptide with given threshold(s)
                if len(all_peptides) > 0:
                    self.peptides[id_].extend({pep: '' for pep in all_peptides
                                            if self._is_good_peptide(pep)})

                miss += 1


    def check_peptide_uniqueness(self):
        """ Check if peptide sequences are unique in the whole fasta file (i.e. 
            has no identical match in other proteins)"""

        # Get all peptides. Use set() to remove peptide duplicates from same protein (if any)
        all_peptides = [peptide for peptides in self.peptides.values()
                        for peptide in set(peptides)]

        # Find frequency of each peptide (does not account peptide duplicates from same protein)
        peptides, count = np.unique(all_peptides, return_counts=True)
        peptides_frequency = dict(zip(peptides, count))

        self.is_unique = {pep: (peptides_frequency[pep] == 1)
                            for peptides in self.peptides.values() for pep in peptides}