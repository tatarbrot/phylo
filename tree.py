#!/usr/bin/python

from Bio import SeqIO
from alignment import *
from msa import *
import numpy as np

class Tree:
    def __init__(self, sequences, names):
        self.sequences = sequences
        self.nodes = list(range(len(self.sequences)))

        self.names = names
        lengths = [len(l) for l in self.sequences]
        lengths = np.diff(lengths)

        self.adjacency_matrix = []

        if any(lengths):
            print('sequences must have same lengths!')
            self.sequences = []


    def infer_tree(self):
        # neighbour joining method
        seq_len = len(self.sequences)

        # start with star topology
        self.nodes.append(len(self.nodes))
        self.clusters = self.nodes[:]

        nodes_len = len(self.nodes)

        d_m = np.zeros([seq_len, seq_len])
        d_m[:,:] = np.inf

        self.adjacency_matrix = np.zeros([nodes_len, nodes_len])

        self.adjacency_matrix[:,:] = np.inf
        self.adjacency_matrix[:,self.adjacency_matrix.shape[1]-1] = 0
        self.adjacency_matrix[self.adjacency_matrix.shape[0]-1, :] = 0
        self.adjacency_matrix[self.adjacency_matrix.shape[0]-1, \
                self.adjacency_matrix.shape[1]-1] = np.inf

        nd = self.nodes[-1]
        child_nodes = np.where(self.adjacency_matrix[nd,:] != np.inf)

        # calculate distances
        for child_node in child_nodes[0]:
            for dist_node in child_nodes[0]:
                if dist_node != child_node:
                    if d_m[child_node,dist_node] == np.inf:
                        d = self.sequence_distance(child_node, dist_node)
                        print('d:',d)
                        d_m[child_node,dist_node] = d
                        d_m[dist_node,child_node] = d
                else:
                    d_m[child_node,dist_node] = 0


        while d_m.shape[0] > 2:
            # calculate q matrix
            qm = np.zeros(d_m.shape)
            qm[:,:] = np.inf
            for i in range(qm.shape[0]):
                for j in range(i):
                    if i != j:
                        qm[i,j] = (seq_len-2)*d_m[i,j] - sum(d_m[i,:]) -sum(d_m[:,j])
                        qm[j,i] = qm[i,j]
                    else:
                        qm[i,j] = np.inf

            # find smallest q value
            minq = np.amin(qm)
            minq = np.where(qm == minq)
            mini = minq[0][0]
            minj = minq[1][0]

            print(d_m)
            print(qm)

            print(self.clusters)
            print(self.clusters[mini], self.clusters[minj])
            print(mini, minj)
            print('--')

            diu = 0.5*d_m[mini,minj] + 1/(2*(seq_len-2))*(sum(d_m[:,minj]) - sum(d_m[mini,:]))
            dju = d_m[mini,minj] - diu

            print(diu)
            print(dju)
            print('...')

            # create new node
            self.nodes.append(len(self.nodes))
            u = self.nodes[-1]

            # update adjacencies
            self.expand_adjacency_matrix()
            print('nd: {} u: {}'.format(nd, u))
            self.adjacency_matrix[u,nd] = 0
            self.adjacency_matrix[nd,u] = 0

            self.adjacency_matrix[u,self.clusters[mini]] = diu
            self.adjacency_matrix[self.clusters[mini],u] = diu
            self.adjacency_matrix[self.clusters[minj],u] = dju
            self.adjacency_matrix[u,self.clusters[minj]] = dju

            # remove old ajdacencies
            self.adjacency_matrix[nd,self.clusters[mini]] = np.inf
            self.adjacency_matrix[self.clusters[mini],nd] = np.inf
            self.adjacency_matrix[nd,self.clusters[minj]] = np.inf
            self.adjacency_matrix[self.clusters[minj],nd] = np.inf

            # start clustering
            if mini < minj:
                n_from = minj
                n_to = mini
            else:
                n_from = mini
                n_to = minj

            self.clusters.pop(n_from)
            self.clusters[n_to] = u

            # update distance matrix
            new_dm = np.delete(d_m, n_from, axis=1)
            new_dm = np.delete(new_dm, n_from, axis=0)
            for i in range(new_dm.shape[0]):
                if n_to != i:
                    if i >= n_from:
                        idx = i+1
                    else:
                        idx = i

                    new_dm[i,n_to] = 0.5*(d_m[idx,mini]+d_m[idx,minj]-d_m[mini,minj])
                    new_dm[n_to, i] = new_dm[i, n_to]
                    print('new_dm[{}, {}] = {}'.format(i,n_to, new_dm[i,n_to]))

            d_m = new_dm

        print(self.adjacency_matrix)

        nw = self.newick_format(nd, [])
        print(nw)

    def newick_format(self, node, ignore):
        ignore.append(node)
        search = self.adjacency_matrix[node, :]
        neighbours = np.where(search != np.inf)
        nw = []
        for nb in neighbours[0]:
            if not nb in ignore:
                if nb < len(self.sequences):
                    print(self.names[nb])
                    nw.append('{}:{}'.format(self.names[nb], search[nb]))
                else:
                    nw.append(self.newick_format(nb, ignore))

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

        return d/l



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

seq_list = []
name_list = []

for i in individuals:
    print(i['name'], i['its'])
    record = SeqIO.read('{}_long.fasta'.format(i['its']), 'fasta')
    seq_list.append(record.seq)
    name_list.append(i['name'])

m = Msa(seq_list, name_list)
sequences, names = m.align()

t = Tree(sequences, names)
t.infer_tree()
