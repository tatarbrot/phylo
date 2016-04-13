#!/usr/bin/env python

import numpy as np

class Alignment:
    def __init__(self, sequences, delta = -2, tp = 'dna'):
        # cost for indel
        self.delta = delta

        self.tp = tp

        self.dna_scores = np.array([\
                ['.', 'a', 't', 'g', 'c'], \
                ['a', 1., -1., -0.5, -1.], \
                ['c', -1., 1., -1., -0.5], \
                ['t', -0.5, -1., 1., -1.], \
                ['g', -1., -0.5, -1., 1.], \
                ])

        self.sequences = ['' for _ in range(len(sequences))]
        for i in range(len(sequences)):
            self.sequences[i] = '-{0}'.format(sequences[i])


        self.alignment_matrix = np.zeros([len(seq) for \
                seq in self.sequences])



    def fill_borders(self):
        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]
        for r in range(rows):
            self.alignment_matrix[r][0] = r*self.delta

        for c in range(cols):
            self.alignment_matrix[0][c] = c*self.delta


    def score(self, row, col):
        r = np.where(self.dna_scores[:][0] == self.sequences[0][row].lower())
        c = np.where(self.dna_scores[0][:] == self.sequences[1][col].lower())

        if len(r[0]) > 0 and len(c[0] > 0):
            return float(self.dna_scores[r[0][0]][c[0][0]])
        else:
            # neutral score
            return 0.0

    def score_cell(self, r, c):
        scores = []
        pos = []
        scores.append(self.alignment_matrix[r-1][c-1] + self.score(r,c))
        pos.append((r-1,c-1))

        scores.append(self.alignment_matrix[r][c-1] + self.delta)
        pos.append((r,c-1))

        scores.append(self.alignment_matrix[r-1][c] + self.delta)
        pos.append((r-1,c))

        # may be biased !
        m = max(scores)
        idx = np.where(scores == m)
        idx =  idx[0][0]

        self.route[r][c] = pos[idx]
        self.alignment_matrix[r][c] = max(scores)

    def align(self):
        self.fill_borders()
        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]
        self.route = np.array([[(0,0) for _ in range(cols)] for _ in range(rows)])
        for r in range(1,rows):
            for c in range(1, cols):
                self.score_cell(r,c)


    def output(self):
        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]

        r = rows-1
        c = cols-1

        s1 = ''
        s2 = ''
        score = 0

        while (r,c) != (0,0):
            p = self.route[r][c]
            score += self.alignment_matrix[r][c]
            if p[0] == r-1 and p[1] == c-1:
                s1 = '{0}{1}'.format(self.sequences[0][r], s1)
                s2 = '{0}{1}'.format(self.sequences[1][c], s2)
            elif p[0] == r and p[1] == c-1:
                s1 = '-{0}'.format(s1)
                s2 = '{0}{1}'.format(self.sequences[1][c], s2)
            elif p[0] == r-1 and p[1] == c:
                s1 = '{0}{1}'.format(self.sequences[0][r], s1)
                s2 = '-{0}'.format(s2)

            r = p[0]
            c = p[1]

        print score
        print 'seq1: ', s1
        print 'seq2: ', s2

    def count_score_from(self, r, c):
        if r == 0 and c == 0:
            return [0],[(0,0)]
        else:
            if c > 0 and r > 0:
                pos = [(r-1, c-1), (r-1, c), (r, c-1)]
            elif c == 0:
                pos = [(r-1, c)]
            elif r == 0:
                pos = [(r, c-1)]

            scores = []
            for p in pos:
                scores.append(self.alignment_matrix[p])

            scores = np.array(scores)
            # find routes
            rts = np.where(scores == np.amax(scores))
            new_scores = []
            new_routes = []
            for i in range(len(rts[0])):
                route = rts[0][i]
                if pos[route] != (0,0):
                    sc, rt = self.count_score_from(pos[route][0], pos[route][1])
                    np_sc = np.array(sc)
                    np_rt = np.array(rt)

                    best_sc = np.where(np_sc == np.amax(sc))
                    for i in range(len(best_sc[0])):
                        new_scores.append(sc[best_sc[0][i]] + scores[route])
                        new_routes.append(rt[best_sc[0][i]])
                        new_routes[-1].append(pos[route])
                else:
                    new_scores.append(scores[route])
                    new_routes.append([pos[route]])


            return new_scores, new_routes



# Tarsius syrichta
s1 = 'AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTATTTTGCCTAGCAAATACAAACTACGAACGAGTCCACAGTCGAACAATAGCACTAGCCCGTGGCCTTCAAACCCTATTACCTCTTGCAGCAACATGATGACTCCTCGCCAGCTTAACCAACCTGGCCCTTCCCCCAACAATTAATTTAATCGGTGAACTGTCCGTAATAATAGCAGCATTTTCATGGTCACACCTAACTATTATCTTAGTAGGCCTTAACACCCTTATCACCGCCCTATATTCCCTATATATACTAATCATAACTCAACGAGGAAAATACACATATCATATCAACAATATCATGCCCCCTTTCACCCGAGAAAATACATTAATAATCATACACCTATTTCCCTTAATCCTACTATCTACCAACCCCAAAGTAATTATAGGAACCATGTACTGTAAATATAGTTTAAACAAAACATTAGATTGTGAGTCTAATAATAGAAGCCCAAAGATTTCTTATTTACCAAGAAAGTA-TGCAAGAACTGCTAACTCATGCCTCCATATATAACAATGTGGCTTTCTT-ACTTTTAAAGGATAGAAGTAATCCATCGGTCTTAGGAACCGAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATTTATTTTCATCCTCCATTTTACTATCACTTACACTCTTAATTACCCCATTTATTATTACAACAACTAAAAAATATGAAACACATGCATACCCTTACTACGTAAAAAACTCTATCGCCTGCGCATTTATAACAAGCCTAGTCCCAATGCTCATATTTCTATACACAAATCAAGAAATAATCATTTCCAACTGACATTGAATAACGATTCATACTATCAAATTATGCCTAAGCTT'
# Lemur catta
s2 = 'AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTATTCTGTCTAGCCAACTCTAACTACGAACGAATCCATAGCCGTACAATACTACTAGCACGAGGGATCCAAACCATTCTCCCTCTTATAGCCACCTGATGACTACTCGCCAGCCTAACTAACCTAGCCCTACCCACCTCTATCAATTTAATTGGCGAACTATTCGTCACTATAGCATCCTTCTCATGATCAAACATTACAATTATCTTAATAGGCTTAAATATGCTCATCACCGCTCTCTATTCCCTCTATATATTAACTACTACACAACGAGGAAAACTCACATATCATTCGCACAACCTAAACCCATCCTTTACACGAGAAAACACCCTTATATCCATACACATACTCCCCCTTCTCCTATTTACCTTAAACCCCAAAATTATTCTAGGACCCACGTACTGTAAATATAGTTTAAA-AAAACACTAGATTGTGAATCCAGAAATAGAAGCTCAAAC-CTTCTTATTTACCGAGAAAGTAATGTATGAACTGCTAACTCTGCACTCCGTATATAAAAATACGGCTATCTCAACTTTTAAAGGATAGAAGTAATCCATTGGCCTTAGGAGCCAAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATCTATTATCCTCTTTCACCCTTGTCACACTGATTATCCTAACTTTACCTATCATTATAAACGTTACAAACATATACAAAAACTACCCCTATGCACCATACGTAAAATCTTCTATTGCATGTGCCTTCATCACTAGCCTCATCCCAACTATATTATTTATCTCCTCAGGACAAGAAACAATCATTTCCAACTGACATTGAATAACAATCCAAACCCTAAAACTATCTATTAGCTT'

# Tarsius syrichta
s1 = 'AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTATTTTGCCTAGCAAATACAAACTACGAACGAGTCCACAGTCGAACAATAGCACTAGCCCGTGGCCTTCAAACCCTATTACCTCTTGCAGCAACATGATGACTCCTCGCCAGCTTAACCAACCTGGCCCTTCCCCCAACAATTAATTTAATCGGTGAACTGTCCGTAATAATAGCAGCATTTTCATGGTCACACCTAACTATTATCTTAGTAGGCCTTAACACCCTTATCACCGCCCTATATTCCCTATATATACTAATCATAACTCAACGAGGAAAATACACATATCATATCAACAATATCATGCCCCCTTTCACCCGAGAAAATACATTAATAATCATACACCTATTTCCCTTAATCCTACTATCTACCAACCCCAAAGTAATTATAGGAACCATGTACTGTAAATATAGTTTAAACAAAACATTAGATTGTGAG'
# Lemur catta
s2 = 'AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTATTCTGTCTAGCCAACTCTAACTACGAACGAATCCATAGCCGTACAATACTACTAGCACGAGGGATCCAAACCATTCTCCCTCTTATAGCCACCTGATGACTACTCGCCAGCCTAACTAACCTAGCCCTACCCACCTCTATCAATTTAATTGGCGAACTATTCGTCACTATAGCATCCTTCTCATGATCAAACATTACAATTATCTTAATAGGCTTAAATATGCTCATCACCGCTCTCTATTCCCTCTATATATTAACTACTACACAACGAGGAAAACTCACATATCATTCGCACAACCTAAACCCATCCTTTACACGAGAAAACACCCTTATATCCATACACATACTCCCCCTTCTCCTATTTACCTTAAACCCCAAAATTATTCTAGGACCCACGTACTGTAAATATAGTTTAAA-AAAACACTAGATTGTGAA'

a = Alignment([s1, s2])
#a = Alignment(['VAHVDDMPNALSALSDLHAHKL', 'AIQLQVTGVVVTDATLKNLGSVHVSKG'])
#a = Alignment(['AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTA', 'AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTA'])
#a = Alignment(['what', 'why'])
#a = Alignment(['GATTAG', 'ATTAC'])
a.align()
a.output()

