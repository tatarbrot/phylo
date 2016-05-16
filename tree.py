#!/usr/bin/python

from Bio import Entrez, SeqIO
from alignment import *
from msa import *
import numpy as np

class Tree:
    def __init__(self, sequences, names):
        self.sequences = sequences
        self.names = names


# data from appendix S1 from Maguilla et al. 2015
individuals = [
        {'name': 'C. arenaria', 'matK': 'KP980114', 'its': 'AY757404'},
        {'name': 'C. bromoides Willd.', 'matK': 'KP980015', 'its': 'AY757404' },
        {'name': 'C. deweyana Schwein.', 'matK': 'KP980016', 'its': 'AY757412'},
        {'name': 'C. leptopoda Mack.', 'matK': 'KP980068', 'its': 'KP980439'},
        {'name': 'C. disperma Dewey', 'matK': 'KP980093', 'its': 'KP980468'},
        {'name': 'C. elongata', 'matK': 'KP980000', 'its': 'KP980370'},
        {'name': 'C. arcta', 'matK': 'KP980047', 'its': 'KP980417'},
        {'name': 'C. arctiformis', 'matK': 'KP980004', 'its': 'KP980375'},
        {'name': 'C. bonanzensis', 'matK': 'KP980038', 'its': 'KP980408'},

        ]

Entrez.email = 'tobias.moser@gmx.ch'
seq_list = []
name_list = []

for i in individuals:
    print(i['name'], i['its'])
    h = Entrez.efetch(db='nucleotide', id=i['its'], \
            rettype='fasta', strand=1, seq_start=0, seq_stop=100)
    record = SeqIO.read(h, 'fasta')
    h.close()
    seq_list.append(record.seq)
    name_list.append(i['name'])

m = Msa(seq_list, name_list)
sequences, names = m.align()

print(sequences, names)
