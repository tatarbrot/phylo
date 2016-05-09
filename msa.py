#!/usr/bin/env python

from alignment import *
import numpy as np
import re

class Msa:
    def __init__(self, sequences, names = []):
        self.sequences = sequences
        self.names = names

        self.nodes = []
        self.clusters = []

        # adjacency matrix
        # indexes with no sequence in self.sequences are
        # nodes inside the tree
        seq_len = len(sequences)
        self.adjacency_matrix = np.zeros([seq_len, seq_len])

        self.root_node = -1

    def align(self):
        # create guid tree
        self.build_guide_tree()

        # root tree
        self.root_tree()

        if self.root_node != -1:
            # continue, finally align sequences
            self.align_sequences(self.root_node, [])
        else:
            print 'could root tree'

    def align_sequences(self, node, ignore):
        # get neighbour_nodes
        search_in

    def get_neighbour_nodes(self, node):
        search_in = self.adjacency_matrix[start][:]
        neighbours = np.where(search_in > 0)

        return neighbours[0]

    def expand_adjacency_matrix(self):
        # expand adjacency matrix
        adjacency_entry = np.zeros([1,self.adjacency_matrix.shape[1]])
        self.adjacency_matrix = np.concatenate((self.adjacency_matrix, adjacency_entry), axis=0)

        adjacency_entry = np.zeros([self.adjacency_matrix.shape[0], 1])
        self.adjacency_matrix = np.concatenate((self.adjacency_matrix, adjacency_entry), axis=1)

    def root_tree(self):
        seq_len = len(self.sequences)
        tree_distances = np.zeros([seq_len, seq_len])
        routes = [[] for _ in range(seq_len)]
        for i in range(seq_len):
            tree_distances[i][:], routes[i] = self.node_distance(i, [], 0, [])

        m = np.amax(tree_distances)
        long_dist = np.where(tree_distances == m)

        d_from = long_dist[0][0]
        d_to = long_dist[1][0]
        root_route = []
        i = 0
        while len(root_route) < 1:
            if routes[d_from][i][-1] == d_to:
                root_route = routes[d_from][i]
            else:
                i += 1

        # add root_node
        self.nodes.append(len(self.nodes))
        self.root_node = self.nodes[-1]

        self.expand_adjacency_matrix()

        mid = m/2
        start_node = -1
        stop_node = -1
        idx = 0
        stop_dist = 0
        start_dist = 0
        while start_node < 0 and idx < len(root_route):
            start = root_route[idx]
            stop = root_route[idx+1]
            stop_dist += self.adjacency_matrix[start][stop]
            if stop_dist > mid:
                start_node = start
                stop_node = stop
            else:
                start_dist = stop_dist

            idx += 1

        # remove distances
        self.adjacency_matrix[start_node][stop_node] = 0
        self.adjacency_matrix[stop_node][start_node] = 0

        # set distances to root_node
        stop_dist = abs(mid-stop_dist)
        start_dist = abs(mid-start_dist)

        self.adjacency_matrix[stop_node][self.root_node] = stop_dist
        self.adjacency_matrix[self.root_node][stop_node] = stop_dist

        self.adjacency_matrix[start_node][self.root_node] = start_dist
        self.adjacency_matrix[self.root_node][start_node] = start_dist

    def node_distance(self, start, ignore, add, route):
        seq_len = len(self.sequences)
        distances = np.zeros([seq_len])

        nb = self.get_neighbour_nodes(start)

        ignore.append(start)
        for n in nb:
            if not n in ignore:
                if n < len(self.sequences):
                    distances[n] += add + self.adjacency_matrix[start][n]
                    rt = ignore[:]
                    rt.append(n)
                    route.append(rt)
                else:
                    a = add + self.adjacency_matrix[start][n]
                    d, route = self.node_distance(n, ignore, a, route)
                    distances = np.add(d, distances)

        return distances, route

    def build_guide_tree(self):

        seq_len = len(self.sequences)

        self.nodes = range(seq_len)
        self.clusters = range(seq_len)

        d_m = np.zeros([seq_len, seq_len])
        while len(self.clusters) > 1:
            for i in range(d_m.shape[0]):
                for j in range(i):
                    if d_m[i][j] == 0:
                        d_m[i][j] = self.average_sequence_distance(i,j)


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

            d_m = np.delete(d_m, (c_from), axis = 0)
            d_m = np.delete(d_m, (c_from), axis = 1)
            d_m[c_to] = 0


    def hamming_distance(self, s0, s1):
        d = 0
        l = len(s0)
        for i in range(l):
            if s0[i] != s1[i]:
                d += 1

        return float(d)/l

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

    def average_sequence_distance(self, i, j):
        # group average & hamming distance

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
                d += self.hamming_distance(aligned_sequences[0][0], \
                        aligned_sequences[1][0])

        return d/(ni*nj)

# Tarsius syrichta
tarsius = 'AAGTTTCATTGGAGCCACCACTCTTATAATTGCCCATGGCCTCACCTCCTCCCTATTATTTTGCCTAGCAAATACAAACTACGAACGAGTCCACAGTCGAACAATAGCACTAGCCCGTGGCCTTCAAACCCTATTACCTCTTGCAGCAACATGATGACTCCTCGCCAGCTTAACCAACCTGGCCCTTCCCCCAACAATTAATTTAATCGGTGAACTGTCCGTAATAATAGCAGCATTTTCATGGTCACACCTAACTATTATCTTAGTAGGCCTTAACACCCTTATCACCGCCCTATATTCCCTATATATACTAATCATAACTCAACGAGGAAAATACACATATCATATCAACAATATCATGCCCCCTTTCACCCGAGAAAATACATTAATAATCATACACCTATTTCCCTTAATCCTACTATCTACCAACCCCAAAGTAATTATAGGAACCATGTACTGTAAATATAGTTTAAACAAAACATTAGATTGTGAGTCTAATAATAGAAGCCCAAAGATTTCTTATTTACCAAGAAAGTA-TGCAAGAACTGCTAACTCATGCCTCCATATATAACAATGTGGCTTTCTT-ACTTTTAAAGGATAGAAGTAATCCATCGGTCTTAGGAACCGAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATTTATTTTCATCCTCCATTTTACTATCACTTACACTCTTAATTACCCCATTTATTATTACAACAACTAAAAAATATGAAACACATGCATACCCTTACTACGTAAAAAACTCTATCGCCTGCGCATTTATAACAAGCCTAGTCCCAATGCTCATATTTCTATACACAAATCAAGAAATAATCATTTCCAACTGACATTGAATAACGATTCATACTATCAAATTATGCCTAAGCTT'
# Lemur catta
lemur = 'AAGCTTCATAGGAGCAACCATTCTAATAATCGCACATGGCCTTACATCATCCATATTATTCTGTCTAGCCAACTCTAACTACGAACGAATCCATAGCCGTACAATACTACTAGCACGAGGGATCCAAACCATTCTCCCTCTTATAGCCACCTGATGACTACTCGCCAGCCTAACTAACCTAGCCCTACCCACCTCTATCAATTTAATTGGCGAACTATTCGTCACTATAGCATCCTTCTCATGATCAAACATTACAATTATCTTAATAGGCTTAAATATGCTCATCACCGCTCTCTATTCCCTCTATATATTAACTACTACACAACGAGGAAAACTCACATATCATTCGCACAACCTAAACCCATCCTTTACACGAGAAAACACCCTTATATCCATACACATACTCCCCCTTCTCCTATTTACCTTAAACCCCAAAATTATTCTAGGACCCACGTACTGTAAATATAGTTTAAA-AAAACACTAGATTGTGAATCCAGAAATAGAAGCTCAAAC-CTTCTTATTTACCGAGAAAGTAATGTATGAACTGCTAACTCTGCACTCCGTATATAAAAATACGGCTATCTCAACTTTTAAAGGATAGAAGTAATCCATTGGCCTTAGGAGCCAAAAA-ATTGGTGCAACTCCAAATAAAAGTAATAAATCTATTATCCTCTTTCACCCTTGTCACACTGATTATCCTAACTTTACCTATCATTATAAACGTTACAAACATATACAAAAACTACCCCTATGCACCATACGTAAAATCTTCTATTGCATGTGCCTTCATCACTAGCCTCATCCCAACTATATTATTTATCTCCTCAGGACAAGAAACAATCATTTCCAACTGACATTGAATAACAATCCAAACCCTAAAACTATCTATTAGCTT'
# homo sapiens
homo = 'AAGCTTCACCGGCGCAGTCATTCTCATAATCGCCCACGGGCTTACATCCTCATTACTATTCTGCCTAGCAAACTCAAACTACGAACGCACTCACAGTCGCATCATAATCCTCTCTCAAGGACTTCAAACTCTACTCCCACTAATAGCTTTTTGATGACTTCTAGCAAGCCTCGCTAACCTCGCCTTACCCCCCACTATTAACCTACTGGGAGAACTCTCTGTGCTAGTAACCACGTTCTCCTGATCAAATATCACTCTCCTACTTACAGGACTCAACATACTAGTCACAGCCCTATACTCCCTCTACATATTTACCACAACACAATGGGGCTCACTCACCCACCACATTAACAACATAAAACCCTCATTCACACGAGAAAACACCCTCATGTTCATACACCTATCCCCCATTCTCCTCCTATCCCTCAACCCCGACATCATTACCGGGTTTTCCTCTTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTTA-CGACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCATGCCCCCATGTCTAACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGCACACTACTATAACCACCCTAACCCTGACTTCCCTAATTCCCCCCATCCTTACCACCCTCGTTAACCCTAACAAAAAAAACTCATACCCCCATTATGTAAAATCCATTGTCGCATCCACCTTTATTATCAGTCTCTTCCCCACAACAATATTCATGTGCCTAGACCAAGAAGTTATTATCTCGAACTGACACTGAGCCACAACCCAAACAACCCAGCTCTCCCTAAGCTT'
# pan
pan = 'AAGCTTCACCGGCGCAATTATCCTCATAATCGCCCACGGACTTACATCCTCATTATTATTCTGCCTAGCAAACTCAAATTATGAACGCACCCACAGTCGCATCATAATTCTCTCCCAAGGACTTCAAACTCTACTCCCACTAATAGCCTTTTGATGACTCCTAGCAAGCCTCGCTAACCTCGCCCTACCCCCTACCATTAATCTCCTAGGGGAACTCTCCGTGCTAGTAACCTCATTCTCCTGATCAAATACCACTCTCCTACTCACAGGATTCAACATACTAATCACAGCCCTGTACTCCCTCTACATGTTTACCACAACACAATGAGGCTCACTCACCCACCACATTAATAACATAAAGCCCTCATTCACACGAGAAAATACTCTCATATTTTTACACCTATCCCCCATCCTCCTTCTATCCCTCAATCCTGATATCATCACTGGATTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGACAACAGAGGCTCA-CGACCCCTTATTTACCGAGAAAGCT-TATAAGAACTGCTAATTCATATCCCCATGCCTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCCATCCGTTGGTCTTAGGCCCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACCATGTATACTACCATAACCACCTTAACCCTAACTCCCTTAATTCTCCCCATCCTCACCACCCTCATTAACCCTAACAAAAAAAACTCATATCCCCATTATGTGAAATCCATTATCGCGTCCACCTTTATCATTAGCCTTTTCCCCACAACAATATTCATATGCCTAGACCAAGAAGCTATTATCTCAAACTGGCACTGAGCAACAACCCAAACAACCCAGCTCTCCCTAAGCTT'
# gorilla
gorilla = 'AAGCTTCACCGGCGCAGTTGTTCTTATAATTGCCCACGGACTTACATCATCATTATTATTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATTCTCTCTCAAGGACTCCAAACCCTACTCCCACTAATAGCCCTTTGATGACTTCTGGCAAGCCTCGCCAACCTCGCCTTACCCCCCACCATTAACCTACTAGGAGAGCTCTCCGTACTAGTAACCACATTCTCCTGATCAAACACCACCCTTTTACTTACAGGATCTAACATACTAATTACAGCCCTGTACTCCCTTTATATATTTACCACAACACAATGAGGCCCACTCACACACCACATCACCAACATAAAACCCTCATTTACACGAGAAAACATCCTCATATTCATGCACCTATCCCCCATCCTCCTCCTATCCCTCAACCCCGATATTATCACCGGGTTCACCTCCTGTAAATATAGTTTAACCAAAACATCAGATTGTGAATCTGATAACAGAGGCTCA-CAACCCCTTATTTACCGAGAAAGCT-CGTAAGAGCTGCTAACTCATACCCCCGTGCTTGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAACTATGTACGCTACCATAACCACCTTAGCCCTAACTTCCTTAATTCCCCCTATCCTTACCACCTTCATCAATCCTAACAAAAAAAGCTCATACCCCCATTACGTAAAATCTATCGTCGCATCCACCTTTATCATCAGCCTCTTCCCCACAACAATATTTCTATGCCTAGACCAAGAAGCTATTATCTCAAGCTGACACTGAGCAACAACCCAAACAATTCAACTCTCCCTAAGCTT'
# pongo
pongo = 'AAGCTTCACCGGCGCAACCACCCTCATGATTGCCCATGGACTCACATCCTCCCTACTGTTCTGCCTAGCAAACTCAAACTACGAACGAACCCACAGCCGCATCATAATCCTCTCTCAAGGCCTTCAAACTCTACTCCCCCTAATAGCCCTCTGATGACTTCTAGCAAGCCTCACTAACCTTGCCCTACCACCCACCATCAACCTTCTAGGAGAACTCTCCGTACTAATAGCCATATTCTCTTGATCTAACATCACCATCCTACTAACAGGACTCAACATACTAATCACAACCCTATACTCTCTCTATATATTCACCACAACACAACGAGGTACACCCACACACCACATCAACAACATAAAACCTTCTTTCACACGCGAAAATACCCTCATGCTCATACACCTATCCCCCATCCTCCTCTTATCCCTCAACCCCAGCATCATCGCTGGGTTCGCCTACTGTAAATATAGTTTAACCAAAACATTAGATTGTGAATCTAATAATAGGGCCCCA-CAACCCCTTATTTACCGAGAAAGCT-CACAAGAACTGCTAACTCTCACT-CCATGTGTGACAACATGGCTTTCTCAGCTTTTAAAGGATAACAGCTATCCCTTGGTCTTAGGATCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAACAGCCATGTTTACCACCATAACTGCCCTCACCTTAACTTCCCTAATCCCCCCCATTACCGCTACCCTCATTAACCCCAACAAAAAAAACCCATACCCCCACTATGTAAAAACGGCCATCGCATCCGCCTTTACTATCAGCCTTATCCCAACAACAATATTTATCTGCCTAGGACAAGAAACCATCGTCACAAACTGATGCTGAACAACCACCCAGACACTACAACTCTCACTAAGCTT'
# hylobates
hylobates = 'AGCTTTACAGGTGCAACCGTCCTCATAATCGCCCACGGACTAACCTCTTCCCTGCTATTCTGCCTTGCAAACTCAAACTACGAACGAACTCACAGCCGCATCATAATCCTATCTCGAGGGCTCCAAGCCTTACTCCCACTGATAGCCTTCTGATGACTCGCAGCAAGCCTCGCTAACCTCGCCCTACCCCCCACTATTAACCTCCTAGGTGAACTCTTCGTACTAATGGCCTCCTTCTCCTGGGCAAACACTACTATTACACTCACCGGGCTCAACGTACTAATCACGGCCCTATACTCCCTTTACATATTTATCATAACACAACGAGGCACACTTACACACCACATTAAAAACATAAAACCCTCACTCACACGAGAAAACATATTAATACTTATGCACCTCTTCCCCCTCCTCCTCCTAACCCTCAACCCTAACATCATTACTGGCTTTACTCCCTGTAAACATAGTTTAATCAAAACATTAGATTGTGAATCTAACAATAGAGGCTCG-AAACCTCTTGCTTACCGAGAAAGCC-CACAAGAACTGCTAACTCACTATCCCATGTATGACAACATGGCTTTCTCAACTTTTAAAGGATAACAGCTATCCATTGGTCTTAGGACCCAAAAATTTTGGTGCAACTCCAAATAAAAGTAATAGCAATGTACACCACCATAGCCATTCTAACGCTAACCTCCCTAATTCCCCCCATTACAGCCACCCTTATTAACCCCAATAAAAAGAACTTATACCCGCACTACGTAAAAATGACCATTGCCTCTACCTTTATAATCAGCCTATTTCCCACAATAATATTCATGTGCACAGACCAAGAAACCATTATTTCAAACTGACACTGAACTGCAACCCAAACGCTAGAACTCTCCCTAAGCTT'

sequences = [tarsius, lemur, homo, pan, gorilla, pongo, hylobates]
names = ['Tarsius syrichta','Lemur catta','Homo sapiens', \
        'Pan','Gorilla','Pongo','Hylobates']
thermotoga = 'ATGGCGAAGGAAAAATTTGTGAGAACAAAACCGCATGTTAACGTTGGAAC'
thermophi = 'ATGGCGAAGGGCGAGTTTGTTCGGACGAAGCCTCACGTGAACGTGGGGAC'
aquaticus = 'ATGGCGAAGGGCGAGTTTATCCGGACGAAGCCCCACGTGAACGTGGGGAC'
deinonema = 'ATGGCTAAGGGAACGTTTGAACGCACCAAACCCCACGTGAACGTGGGCAC'
chlamydia = 'ATGTCAAAAGAAACTTTTCAACGTAATAAGCCTCATATCAACATAGGGGC'
flexistips = 'ATGTCCAAGCAAAAGTACGAAAGGAAGAAACCTCACGTAAACGTAGGCAC'
borrelia = 'ATGGCAAAAGAAGTTTTTCAAAGAACAAAGCCGCACATGAATGTTGGAAC'
sequences = [thermotoga, thermophi, aquaticus, deinonema, chlamydia, \
        flexistips, borrelia]
names = ['thermotoga', 'thermophi', 'aquaticus', 'deinonema', 'chlamydia', \
        'flexistips', 'borrelia']

m = Msa(sequences, names)
m.align()
