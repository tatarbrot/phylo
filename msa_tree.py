#!/usr/bin/env python

from alignment import *
from tree import *
import numpy as np

class Msa(Tree):
    def __init__(self, sequences, names = []):
        Tree.__init__(self, sequences, names)


    def align(self):
        print("Performing Multiptle sequence alignment")
        print("Inferring guide tree")
        # create guid tree
        self.infer_tree()

        # root tree
        self.root_tree()

        print("Aligning sequences")
        if self.root_node != -1:
            # continue, finally align sequences
            seq = self.align_sequences(self.root_node, [])

            sequences = []
            names = []
            for s in seq:
                sequences.append(s[0])
                names.append(self.names[s[1]])
            return sequences, names
        else:
            print('could not root tree')

    def align_sequences(self, node, ignore):
        # get neighbour_nodes
        neighbours = self.get_neighbour_nodes(node)
        ignore.append(node)
        all_seq_len = len(self.sequences)
        tree_nodes = [n < all_seq_len for n in neighbours]

        align_seq = []
        for n in neighbours:
            if not n in ignore:
                if n < all_seq_len:
                    align_seq.append([(self.sequences[n], n)])
                else:
                    align_seq.append(self.align_sequences(n, ignore))

        seq = []
        if len(align_seq) > 2:
            print('unexpected number (>2) of sequences to align')
        elif len(align_seq) == 1:
            return align_seq[0]
        else:
            aseq = [[s[0] for s in c] for c in align_seq]
            a = Alignment(aseq)
            a.align()
            aligned_sequences = a.output()
            for i in range(len(aligned_sequences)):
                for j in range(len(aligned_sequences[i])):
                    seq.append((aligned_sequences[i][j], align_seq[i][j][1]))
            return seq


    def infer_tree(self):
        seq_len = len(self.sequences)

        d_m = np.zeros([seq_len, seq_len])
        d_m[:][:] = np.inf

        while len(self.clusters) > 1:
            for i in range(d_m.shape[0]):
                for j in range(i):
                    if d_m[i,j] == np.inf:
                        d_m[i,j] = self.sequence_distance(i,j)
                        d_m[j,i] = d_m[i,j]


            # find minimums
            minimas = []
            for i in range(1,d_m.shape[0]):
                minimas.append(np.amin(d_m[i][0:i]))

            min_dist = min(minimas)
            min_idx = np.where(d_m == min_dist)

            if min_idx[0][0] < min_idx[1][0]:
                c_from = min_idx[1][0]
                c_to = min_idx[0][0]
            else:
                c_from = min_idx[0][0]
                c_to = min_idx[1][0]

            # create new node
            self.nodes.append(len(self.nodes))
            self.expand_adjacency_matrix()

            # set adjacencies
            self.adjacency_matrix[-1][self.clusters[c_from]] = min_dist/2
            self.adjacency_matrix[-1][self.clusters[c_to]] = min_dist/2

            self.adjacency_matrix[self.clusters[c_from]][-1] = min_dist/2
            self.adjacency_matrix[self.clusters[c_to]][-1] = min_dist/2

            self.clusters[c_to] = self.nodes[-1]
            self.clusters.pop(c_from)

            new_d_m = np.delete(d_m, (c_from), axis = 0)
            new_d_m = np.delete(new_d_m, (c_from), axis = 1)

            for i in range(new_d_m.shape[0]):
                    if i != c_to:
                        if i >= c_from:
                            idx = i+1
                        else:
                            idx = i

                        new_d = (d_m[idx, c_to] + d_m[idx, c_from])/2
                        new_d_m[c_to, i] = new_d
                        new_d_m[i, c_to] = new_d

            d_m = new_d_m

    def tree_members(self, ti, ignore):
        if ti >= len(self.sequences):
            # get neighbours
            mem = self.get_neighbour_nodes(ti)
            members = []
            ignore.append(ti)
            for m in mem:
                if not m in ignore:
                    if m < len(self.sequences):
                        members.append(m)
                    else:
                        members.extend(self.tree_members(m, ignore))
        else:
            members = [ti]

        return members

    def sequence_distance(self, i, j):
        # group average & jaccard distance
        ti = self.clusters[i]
        tj = self.clusters[j]

        mi = self.tree_members(ti, [])
        mj = self.tree_members(tj, [])
        ni = len(mi)
        nj = len(mj)

        d = 0
        for k in range(ni):
            s1 = self.sequences[mi[k]]
            for l in range(nj):
                s2 = self.sequences[mj[l]]
                a = Alignment([[s1], [s2]])
                a.align()
                aligned_sequences = a.output()
                d += self.jaccard_distance(aligned_sequences[0][0], \
                        aligned_sequences[1][0])

        return d/(ni*nj)
