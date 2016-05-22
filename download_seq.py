#!/usr/bin/python
from Bio import Entrez, SeqIO
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
            rettype='fasta', strand=1)
    record = SeqIO.read(h, 'fasta')
    SeqIO.write(record, '{}_long.fasta'.format(i['its']), 'fasta')
    h.close()
