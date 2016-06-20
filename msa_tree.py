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

        # create guide tree
        self.infer_tree()

        # root tree
        self.root_tree()

        print("Aligning sequences")
        # continue if rooting worked
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
        # align sequences according to the inferred guide tree
        # recursive function

        # find all neighbors and add self to ignore list
        neighbours = self.get_neighbour_nodes(node)
        ignore.append(node)
        all_seq_len = len(self.sequences)
        tree_nodes = [n < all_seq_len for n in neighbours]

        # create list with what has to be aligned
        align_seq = []
        for n in neighbours:
            if not n in ignore:
                if n < all_seq_len:
                    # if leaf node, add sequences to aligning list
                    align_seq.append([(self.sequences[n], n)])
                else:
                    # if internal node, first align the children of the
                    # according node
                    align_seq.append(self.align_sequences(n, ignore))


        # align what has been added to the alignment list
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
        # create guide tree with UPGMA method

        seq_len = len(self.sequences)

        # create distance matrix
        d_m = np.zeros([seq_len, seq_len])
        d_m[:][:] = np.inf

        # continue until there is only one cluster left
        while len(self.clusters) > 1:
            # update distnaces
            for i in range(d_m.shape[0]):
                for j in range(i):
                    if d_m[i,j] == np.inf:
                        d_m[i,j] = self.sequence_distance(i,j)
                        d_m[j,i] = d_m[i,j]


            # find minima
            minimas = []
            for i in range(1,d_m.shape[0]):
                minimas.append(np.amin(d_m[i][0:i]))

            min_dist = min(minimas)
            min_idx = np.where(d_m == min_dist)

            # find smaller index and join nodes into smaller
            # index
            if min_idx[0][0] < min_idx[1][0]:
                c_from = min_idx[1][0]
                c_to = min_idx[0][0]
            else:
                c_from = min_idx[0][0]
                c_to = min_idx[1][0]

            # create new internal node in adjacency matrix
            self.nodes.append(len(self.nodes))
            self.expand_adjacency_matrix()

            # set distances in adjacency matrix
            self.adjacency_matrix[-1][self.clusters[c_from]] = float(min_dist)/2
            self.adjacency_matrix[-1][self.clusters[c_to]] = float(min_dist)/2
            self.adjacency_matrix[self.clusters[c_from]][-1] = float(min_dist)/2
            self.adjacency_matrix[self.clusters[c_to]][-1] = float(min_dist)/2

            # update cluster with smaller index
            self.clusters[c_to] = self.nodes[-1]
            self.clusters.pop(c_from)

            # update distance matrix
            new_d_m = np.delete(d_m, (c_from), axis = 0)
            new_d_m = np.delete(new_d_m, (c_from), axis = 1)
            for i in range(new_d_m.shape[0]):
                    if i != c_to:
                        if i >= c_from:
                            idx = i+1
                        else:
                            idx = i

                        new_d = float(d_m[idx, c_to] + d_m[idx, c_from])/2
                        new_d_m[c_to, i] = new_d
                        new_d_m[i, c_to] = new_d

            d_m = new_d_m

    def tree_members(self, ti, ignore):
        # find children of current node and return list of
        # leaf nodes
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

        return float(d)/(ni*nj)
