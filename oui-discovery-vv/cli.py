#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
"""
Function used to parse argument when digestool is run with CLI
"""
# ------------------------------------------------------------------------------
import argparse
# ------------------------------------------------------------------------------

def parse_args():
    """ Parse command line arguments to run a cleave protein of the given 
        fasta file"""

    parser = argparse.ArgumentParser(description='Digest proteins of a given fasta file')
    parser.add_argument(
            'fasta_file',
            metavar='</path/to/file.fasta>'
    )
    parser.add_argument(
            '-l', '--min_length',
            metavar='',
            default=7,
            help='minimal length of peptide sequences (default: 7)'
    )
    parser.add_argument(
            '-L', '--max_length',
            metavar='',
            default=None,
            help='maximal length of peptide sequences (default: None)'
    )
    parser.add_argument(
            '-M', '--miscleavages',
            metavar='',
            default=2,
            help='maximum of miscleavages allowed (default: 2)'
    )
    parser.add_argument(
            '-m', '--mass',
            metavar='',
            default=4600,
            help='maximal molecular mass of peptide sequences (default: 4600 dalton)'
    )
    parser.add_argument(
            '-u ', '--unique',
            action='store_true',
            help='add flag to inform if the peptides are unique or not'
    )
    parser.add_argument(
            '-e', '--enzyme',
            metavar='',
            default='trypsin',
            help='enzyme used to perform digestion (default: trypsin)'
    )
    parser.add_argument(
            '-o', '--output',
            metavar='',
            help='output directory (default: /path/to/file.fasta_digestedPeptides.csv)'
    )

    args = parser.parse_args()
    return args
