#!/usr/bin/python

from Bio import Entrez, SeqIO

individuals = [
        {'name': 'C. arenaria', 'matK': 'KP980114', 'its': 'AY757404'},
        {'name': 'C. bromoides Willd.', 'matK': 'KP980015', 'its': 'AY757404' },
        {'name': 'C. deweyana Schwein.', 'matK': 'KP980016', 'its': 'AY757412'},
        {'name': 'C. leptopoda Mack.', 'matK': 'KP980068', 'its': 'KP980439'},
        {'name': 'C. disperma Dewey', 'matK': 'KP980093', 'its': 'KP980468'},
        {'name': 'C. elongata L.', 'matK': 'KP980000', 'its': 'KP980370'},
        {'name': 'C. arcta Boott', 'matK': 'KP980047', 'its': 'KP980417'},
        {'name': 'C. arctiformis Mack', 'matK': 'KP980004', 'its': 'KP980375'},
        {'name': 'C. bonanzensis Britton', 'matK': 'KP980038', 'its': 'KP980408'},
        {'name': 'C. brunnescens [Pers.] Poir.', 'matK': 'KP980130', 'its': 'KP980514'},
        {'name': 'C. canescens L.', 'matK': '', 'its': 'KP980449'},
        {'name': 'C. diastena V.I.Krecz.', 'matK': 'KP980089', 'its': 'KP980464'},
        {'name': 'C. furva Webb', 'matK': 'KP980138', 'its': 'KP980528'},
        {'name': 'C. glareosa Schkuhr ex', 'matK': 'KP980072', 'its': 'KP980444'},
        {'name': 'C. heleonastes Ehrh. ex L.f.', 'matK': 'KP980040', 'its': 'KP980410'},
        {'name': 'C. kreczetoviczii T.V.Egor.', 'matK': 'KP980154', 'its': 'KP980546'},
        {'name': 'C. lachenalii Schkuhr', 'matK': 'KP979994', 'its': 'KP980364'},
        {'name': 'C. laeviculmis Meinsh.', 'matK': 'KP980018', 'its': 'KP980388'},
        {'name': 'C. lapponica O.Lang', 'matK': 'KP980074', 'its': 'KP980446'},
        {'name': 'C. loliacea L.', 'matK': 'KP980118', 'its': 'KP980500'},
        {'name': 'C. mackenziei V.I.Krecz.', 'matK': 'KP980071', 'its': 'KP980443'},
        {'name': 'C. marina Dewey', 'matK': 'KP979995', 'its': 'KP980365'},
        {'name': 'C. nemurensis Franch.', 'matK': 'KP980033', 'its': 'KP980403'},
        {'name': 'C. praeceptorum Mack.', 'matK': 'KP980104', 'its': 'KP980481'},
        {'name': 'C. pseudololiacea F.Schmidt.', 'matK': 'KP980088', 'its': 'KP980463'},
        {'name': 'C. tenuiflora Wahlenb.', 'matK': 'KP980045', 'its': 'KP980415'},
        {'name': 'C. traiziscana F.Schmidt.', 'matK': 'KP980115', 'its': 'KP980497'},
        {'name': 'C. trisperma Dewey', 'matK': 'KP980142', 'its': 'KP980532'},
        {'name': 'C. ursina Dewey', 'matK': 'KP980010', 'its': 'KP980380'},
        {'name': 'C. paniculata L.', 'matK': 'KP980057', 'its': 'KP980427'},
        {'name': 'C. illota L.H.Bailey', 'matK': 'KP980017', 'its': 'KP980387'},
        {'name': 'C. spicata Huds.', 'matK': 'KP980082', 'its': 'KP980454'},
        {'name': 'C. remota L.', 'matK': 'KP980019', 'its': 'KP980389'},
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
