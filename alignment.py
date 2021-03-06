#!/usr/bin/env python

import numpy as np

class Alignment:
    def __init__(self, sequences, weights = [], delta = -2, tp = 'dna'):
        # cost for indel
        self.delta = delta

        if len(weights) < 1:
            # equal weights
            self.weights = [[1 for _ in range(len(sequences[i]))] \
                    for i in range(2)]
        else:
            self.weights = weights

        self.weights = np.array(self.weights)

        # normalize weights
        for i in range(len(self.weights)):
            self.weights[i] = np.true_divide(self.weights[i], \
                    np.sum(self.weights[i]))

        self.sequences = sequences
        self.tp = tp

        if self.tp == 'dna':
            self.scoring_table = np.array([\
                    ['#', 'a', 't', 'g', 'c', '.'], \
                    ['a', 1., -1., -0.5, -1., -1], \
                    ['c', -1., 1., -1., -0.5, -1], \
                    ['t', -0.5, -1., 1., -1., -1], \
                    ['g', -1., -0.5, -1., 1., -1], \
                    ['.', -1, -1, -1, -1, 1], \
                    ])

        self.sequences = [['' for _ in range(len(sequences[i]))] for i in range(2)]
        for i in range(len(sequences)):
            for j in range(len(sequences[i])):
                self.sequences[i][j] = '-{0}'.format(sequences[i][j])


        self.alignment_matrix = np.zeros([len(seq[0]) for \
                seq in self.sequences])



    def fill_borders(self):
        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]
        for r in range(rows):
            self.alignment_matrix[r][0] = r*self.delta

        for c in range(cols):
            self.alignment_matrix[0][c] = c*self.delta


    def score(self, row, col, s0 = 0, s1 = 0):
        r = np.where(self.scoring_table[:][0] == self.sequences[0][s0][row].lower())
        c = np.where(self.scoring_table[0][:] == self.sequences[1][s1][col].lower())

        if len(r[0]) > 0 and len(c[0] > 0):
            return float(self.scoring_table[r[0][0]][c[0][0]])
        else:
            # neutral score
            return 0.0

    def score_cell(self, r, c):
        scores = []
        pos = []

        l0 = len(self.sequences[0])
        l1 = len(self.sequences[1])
        match_score = 0

        # get scores of pairwise character alignments
        match_scores = np.zeros(l1*l0)
        for i in range(l0):
            for j in range(l1):
                match_scores[j+i*l1] = \
                        self.score(r,c,i,j)*self.weights[0][i]*self.weights[1][j]

        # create mean of scores
        match_score = np.mean(match_scores)

        # append all possible scores to an array
        # with the according position of the cell
        scores.append(self.alignment_matrix[r-1][c-1] + match_score)
        pos.append((r-1,c-1))

        scores.append(self.alignment_matrix[r][c-1] + self.delta)
        pos.append((r,c-1))

        scores.append(self.alignment_matrix[r-1][c] + self.delta)
        pos.append((r-1,c))

        # get best score
        pos = np.array(pos)
        scores = np.array(scores)

        # get best score
        # pick the first if there are more than one equal
        # to the best score
        m = np.amax(scores)
        idx = np.where(scores == m)
        if len(idx[0]) == 1:
            p = pos[idx[0][0]]
        elif len(idx[0]) > 1:
            check_pos = pos[idx]
            check_scores = []
            for cp in check_pos:
                check_scores.append(self.alignment_matrix[cp[0]][cp[1]])

            new_idx = np.where(np.array(check_scores) == np.amax(check_scores))
            p = check_pos[new_idx[0][0]]

        # set route for the according cell
        self.route[r][c] = p
        # update score of cell
        self.alignment_matrix[r][c] = max(scores)

    def align(self):
        # fill borders with either always deletion/insertion gap penalty
        self.fill_borders()

        # get shape of alignment matrix
        rows = self.alignment_matrix.shape[0]
        cols = self.alignment_matrix.shape[1]
        self.route = np.array([[(0,0) for _ in range(cols)] for _ in range(rows)])

        # create score for each cell
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

        l0 = len(self.sequences[0])
        l1 = len(self.sequences[1])

        s0 = ['' for _ in range(l0)]
        s1 = ['' for _ in range(l1)]


        # walk along the created route and
        # put together the sequences according
        # to the choices made in score_cell()
        while (r,c) != (0,0):
            p = self.route[r][c]
            score += self.alignment_matrix[r][c]

            if p[0] == r-1 and p[1] == c-1:
                # align characters
                for i in range(l0):
                    s0[i] = '{0}{1}'.format(self.sequences[0][i][r], s0[i])
                for j in range(l1):
                    s1[j] = '{0}{1}'.format(self.sequences[1][j][c], s1[j])
            elif p[0] == r and p[1] == c-1:
                # deletion in sequence group 0
                for i in range(l0):
                    s0[i] = '.{0}'.format(s0[i])
                for j in range(l1):
                    s1[j] = '{0}{1}'.format(self.sequences[1][j][c], s1[j])
            elif p[0] == r-1 and p[1] == c:
                # insertion in sequence group 0
                for i in range(l0):
                    s0[i] = '{0}{1}'.format(self.sequences[0][i][r], s0[i])
                for j in range(l1):
                    s1[j] = '.{0}'.format(s1[j])

            r = p[0]
            c = p[1]

        return [s0, s1]
