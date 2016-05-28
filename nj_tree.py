#!/usr/bin/python

from Bio import SeqIO, Phylo
from io import StringIO
from alignment import *
from msa_tree import *
import sys
import numpy as np

class NjTree(Tree):
    def __init__(self, sequences, names):
        Tree.__init__(self, sequences, names)

        lengths = [len(l) for l in self.sequences]
        lengths = np.diff(lengths)

        if any(lengths):
            print('sequences must have same lengths!')
            self.sequences = []


    def create_tree(self):

        self.infer_tree()
        self.root_tree()
        nw = self.newick_format(self.root_node, [])
        return nw

    def infer_tree(self):

        print('Inferring Neighbour-joining tree...')

        # neighbour joining method
        seq_len = len(self.sequences)

        # start with star topology
        self.nodes.append(len(self.nodes))
        self.clusters = self.nodes[:]

        nodes_len = len(self.nodes)

        d_m = np.zeros([seq_len, seq_len])
        d_m[:,:] = np.inf

        self.expand_adjacency_matrix()
        l_entry = [0,0]
        l_entry[0] = self.adjacency_matrix.shape[0]-1
        l_entry[1] = self.adjacency_matrix.shape[1]-1

        # create star topology
        self.adjacency_matrix[:,:] = np.inf
        self.adjacency_matrix[:,l_entry[1]] = 0
        self.adjacency_matrix[l_entry[0], :] = 0
        self.adjacency_matrix[l_entry[0],l_entry[1]] = np.inf

        # get main node of star topology
        nd = self.nodes[-1]
        child_nodes = np.where(self.adjacency_matrix[nd,:] != np.inf)

        # calculate distances
        for child_node in child_nodes[0]:
            for dist_node in child_nodes[0]:
                if dist_node != child_node:
                    if d_m[child_node,dist_node] == np.inf:
                        d = self.sequence_distance(child_node, dist_node)
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

            diu = 0.5*d_m[mini,minj] + float(1)/(2*(seq_len-2))*(sum(d_m[:,minj]) - sum(d_m[mini,:]))
            dju = d_m[mini,minj] - diu

            # create new node
            self.nodes.append(len(self.nodes))
            u = self.nodes[-1]

            #update adjacencies
            self.expand_adjacency_matrix()
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

            #update distance matrix
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

            d_m = new_dm

# data from appendix S1 from Maguilla et al. 2015
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

seq_list = []
name_list = []

lof = 'short'
if len(sys.argv) > 1:
    lof = sys.argv[1]

for i in individuals:
    print(i['name'], i['its'])
    record = SeqIO.read('{}_{}.fasta'.format(i['its'], lof), 'fasta')
    seq_list.append(record.seq)
    name_list.append(i['name'])

# perform alignment
m = Msa(seq_list, name_list)
sequences, names = m.align()

# Infer tree
t = NjTree(sequences, names)
nw = t.create_tree()

print(nw)
handle = StringIO(nw)
disp_tree = Phylo.read(handle, 'newick')
Phylo.draw_ascii(disp_tree)
