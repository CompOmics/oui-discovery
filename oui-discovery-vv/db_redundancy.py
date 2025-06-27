from itertools import product
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
import copy
import pickle

def kmer_aa(word_size=8):
    """
        lower bond of tryptic peptides length https://www.ptglab.com/news/blog/mass-spec-essentials/
    """
    alphabet="ACDEFGHIKLMNPQRSTVWY"
    kmers_seq = []
    for y in product(alphabet, repeat=word_size):
        kmers_seq.append(''.join(y))
    with open('kmers_seq_empty.pkl', 'wb') as f:
        pickle.dump(kmers_seq, f)
        
def db_redundancy(db_fasta,name_pkl,word_size=8): #kmers_seq_empty
    #kmer_dict={k:0 for k in kmers_seq}
    kmer_dict={}
    for record in db_fasta:
        seq=str(db_fasta[record].seq)
        for i in range(0,len(seq)-word_size):
            kmer=seq[i:i+word_size]
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer]=0
    with open(f'{name_pkl}.pickle', 'wb') as handle:
        pickle.dump(kmer_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def get_fasta(dir="./ionbot_openprot/"):
    UPfasta_search=f'{dir}Human_2023_01_canonical.fasta'
    input_file = open(UPfasta_search)
    UPfasta = SeqIO.to_dict(SeqIO.parse(UPfasta_search, "fasta"))
    
    OPfasta_file=f'{dir}openprot_2_0_0-human-ensembl106+refseq+uniprot2022_06_01.fasta'
    input_file = open(OPfasta_file)
    OPfasta = SeqIO.to_dict(SeqIO.parse(OPfasta_file, "fasta"))
    
    TRfasta_search=f'{dir}Human_2023_01_trembl.fasta'
    input_file = open(TRfasta_search)
    TRfasta = SeqIO.to_dict(SeqIO.parse(TRfasta_search, "fasta"))
    
    return UPfasta, OPfasta, TRfasta