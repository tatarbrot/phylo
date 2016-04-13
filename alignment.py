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
        scores.append(self.alignment_matrix[r-1][c-1] + self.score(r,c))
        scores.append(self.alignment_matrix[r][c-1] + self.delta)
        scores.append(self.alignment_matrix[r-1][c] + self.delta)

        print scores, max(scores)
        self.alignment_matrix[r][c] = max(scores)

    def align(self):
        self.fill_borders()

        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]
        for r in range(1,rows):
            for c in range(1, cols):
                self.score_cell(r,c)

        # find routes
        scores, routes = self.count_score_from(rows-1, cols-1)
        self.scores = scores
        self.routes = routes

    def output(self):
        r_old = -1
        c_old = -1

        for i in range(len(self.scores)):
            s1 = ''
            s2 = ''
            for j in range(len(self.routes[i])-1, 0, -1):
                p = self.routes[i][j]
                r = p[0]
                c = p[1]
                if r == r_old:
                    s1 = '-{0}'.format(s1)
                else:
                    s1 = '{0}{1}'.format(self.sequences[0][r], s1)

                if c == c_old:
                    s2 = '-{0}'.format(s2)
                else:
                    s2 = '{0}{1}'.format(self.sequences[1][c], s2)

                r_old = r
                c_old = c

            print self.scores[i]
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

#a = Alignment([s1, s2])
#a = Alignment(['VAHVDDMPNALSALSDLHAHKL', 'AIQLQVTGVVVTDATLKNLGSVHVSKG'])
#a = Alignment(['AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTA', 'AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTA'])
a = Alignment(['why', 'what'])
a.align()
a.output()

