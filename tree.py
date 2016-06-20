#!/usr/bin/python

from alignment import *
import numpy as np
import re

class Tree:
    def __init__(self, sequences, names):
        # create list of ...
        # seqeunces
        self.sequences = sequences
        # nodes
        self.nodes = list(range(len(self.sequences)))
        # clusters
        self.clusters = self.nodes[:]
        # names
        self.names = names

        # create adjacency matrix
        seq_len = len(sequences)
        self.adjacency_matrix = np.zeros([seq_len, seq_len])

        # set adjacencies to infinity
        self.adjacency_matrix[:][:] = np.inf

        # no root node known
        self.root_node = -1

    def infer_tree(self):
        # must be implemented by children
        pass

    def newick_format(self, node, ignore):
        # create newick output
        ignore.append(node)
        search = self.adjacency_matrix[node, :]
        neighbours = np.where(search != np.inf)
        nw = []
        for nb in neighbours[0]:
            if not nb in ignore:
                if nb < len(self.sequences):
                    nam = re.sub('\s','_',self.names[nb])
                    print(nam)
                    nw.append('{}:{}'.format(nam, search[nb]))
                else:
                    nw.append('{}:{}'.format(self.newick_format(nb, ignore), search[nb]))

        nw_str = '('
        for n in nw:
            nw_str += '{},'.format(n)

        nw_str = nw_str[0:len(nw_str)-1]
        nw_str += ')'
        return nw_str


    def expand_adjacency_matrix(self):
        adj_entry = np.zeros([1,self.adjacency_matrix.shape[1]])
        adj_entry[:] = np.inf
        self.adjacency_matrix = np.concatenate((self.adjacency_matrix, adj_entry), axis=0)

        adj_entry = np.zeros([self.adjacency_matrix.shape[0],1])
        adj_entry[:] = np.inf
        self.adjacency_matrix = np.concatenate((self.adjacency_matrix, adj_entry), axis=1)

    def sequence_distance(self, ni, nj):
        # jaccards distance between two sequences
        s1 = self.sequences[ni]
        s2 = self.sequences[nj]

        a = Alignment([[s1],[s2]])
        a.align()
        aligned_sequences = a.output()
        d = self.jaccard_distance(aligned_sequences[0][0], aligned_sequences[1][0])

        return d

    def jaccard_distance(self, s1, s2):
        d = 0
        l = len(s1)
        for i in range(l):
            if s1[i] != s2[i]:
                d += 1

        return float(d)/l

    def get_neighbour_nodes(self, node):
        # return immediate neighbour nodes (not recursive)
        search_in = self.adjacency_matrix[node][:]
        neighbours = np.where(search_in < np.inf)

        return neighbours[0]

    def node_distance(self, start, ignore, add, route):
        # distance between two nodes
        # walk along adjacencies
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

    def root_tree(self):
        print("Rooting tree (mid-point)")
        seq_len = len(self.sequences)
        tree_distances = np.zeros([seq_len, seq_len])
        routes = [[] for _ in range(seq_len)]
        # calculate distances of leaf nodes from the tree distances
        for i in range(seq_len):
            tree_distances[i][:], routes[i] = self.node_distance(i, [], 0, [])

        # take maximal length
        m = np.amax(tree_distances)
        long_dist = np.where(tree_distances == m)

        d_from = long_dist[0][0]
        d_to = long_dist[1][0]

        root_route = []
        i = 0
        # find route which leads from d_from to d_to
        while len(root_route) < 1:
            if routes[d_from][i][-1] == d_to:
                root_route = routes[d_from][i]
            else:
                i += 1

        # add root_node in the middle in between
        self.nodes.append(len(self.nodes))
        self.root_node = self.nodes[-1]
        self.expand_adjacency_matrix()
        # find middle
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
        self.adjacency_matrix[start_node][stop_node] = np.inf
        self.adjacency_matrix[stop_node][start_node] = np.inf

        # set distances to root_node
        stop_dist = abs(mid-stop_dist)
        start_dist = abs(mid-start_dist)
        self.adjacency_matrix[stop_node][self.root_node] = stop_dist
        self.adjacency_matrix[self.root_node][stop_node] = stop_dist
        self.adjacency_matrix[start_node][self.root_node] = start_dist
        self.adjacency_matrix[self.root_node][start_node] = start_dist
