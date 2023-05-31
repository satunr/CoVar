import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt
import collections
import stat
import pickle
import math
import pandas as pd
# import mygene

from networkx.algorithms.community.centrality import girvan_newman
from Nearest import *
from variation import *
from scipy import stats
from copy import deepcopy
from itertools import permutations, combinations
from Genie import GENIE3, create_graph


def maps(g_names):
    M = {}
    L = ['symbol', 'name']
    mg = mygene.MyGeneInfo()

    for gene in g_names:
        try:
            g = mg.getgene(gene)
            g = {field: g[field] for field in g.keys() if field in L}
            M[gene] = g['symbol']
        except:
            M[gene] = gene
        print (gene, M[gene])

    return M


def find_cores(G, cut_off = 5):
    G.remove_edges_from([(u, v) for (u, v) in G.edges() if G[u][v]['weight'] < cut_off])

    I = list(nx.isolates(G))
    G.remove_nodes_from(I)

    labels, lC, G = comm2(G, None)
    Cores = [[gene for gene in labels.keys() if labels[gene] == i] for i in range(lC)]

    return Cores, labels


def comm2(G, g_names):
    H = G.to_undirected()
    comp = girvan_newman(H)
    C = [sorted(c) for c in next(comp)]
    c = 0
    labels = {}

    for l in C:
        for gene in l:
            labels[gene] = c
        c = c + 1

    Z = nx.DiGraph()
    for gene in G.nodes():
        Z.add_node(gene, weight = labels[gene])
    Z.add_edges_from(list(G.edges()))

    return labels, c, Z


how_many_runs = 12
# mg = mygene.MyGeneInfo()

# Combined analysis
CombinedKnn = nx.DiGraph()
Maps = {}
C = []
for ind in range(how_many_runs):
    print ('Iteration 0',
           ind)

    # < Control, Disease, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >
    B = pickle.load(open('GList-Macro-a-' + str(ind) + '.p', 'rb'))
    [_, _, gnames, _, Va, Knn, _, C] = deepcopy(B)
    print ('Iteration 1', ind)

    # M = maps(gnames)

    Cores = []
    for gene_list in C:
        print (gene_list)
        for gene in gene_list:
            Cores.append(gene)

    for u in Knn.nodes():
        if u not in C:
            C.append(u)

    for e in Knn.edges():
        if e not in CombinedKnn.edges():
            CombinedKnn.add_edge(e[0], e[1], weight = 1.0)
        else:
            CombinedKnn[e[0]][e[1]]['weight'] += 1.0

    for gene in Knn.nodes():
        if gene not in Maps.keys():
            # KNN
            Maps[gene] = [0, 1, 0]
        else:
            Maps[gene] = [Maps[gene][0], Maps[gene][1] + 1, Maps[gene][2]]

        if gene in Cores:
            # Cores
            Maps[gene] = [Maps[gene][0], Maps[gene][1], Maps[gene][2] + 1]

        if gene in Va:
            # Variational
            Maps[gene] = [Maps[gene][0] + 1, Maps[gene][1], Maps[gene][2]]

    nx.write_gml(CombinedKnn, 'Comb_IBD.gml')
    pickle.dump(Maps, open('Maps.p', 'wb'))
exit(1)


# Remap nodes (Ensemble to symbol)
# --------------------------------
'''
M = maps(gnames)
pickle.dump(M, open('ens_symbol.p', 'wb'))
'''
'''
# Create Combined network
# -----------------------
Maps = pickle.load(open('Maps.p', 'rb'))
gnames = list(Maps.keys())

M = pickle.load(open('ens_symbol.p', 'rb'))

G = nx.read_gml('Comb_IBD.gml')
G = nx.relabel_nodes(G, {gene: M[gene] for gene in gnames})
print (list(G.nodes()))
nx.write_gml(G, 'Comb_IBD_renamed.gml')

# Find cores in combined network
Cores, labels = find_cores(G)

condense = []
for each in Cores:
    for gene in each:
        condense.append(gene)

print ([len(c) for c in Cores])
'''
'''
# Create spreadsheet
# ------------------
D = {'Ensemble': [], 'Symbol': [], 'KNN': [], 'Variational': [],
     'Core': [], 'Final Core': [], 'Cluster ID': []}

for gene in Maps.keys():
    if M[gene] in labels.keys():
        l = labels[M[gene]]
    else:
        l = -1

    if M[gene] in condense:
        c = 1
    else:
        c = 0

    D['Ensemble'].append(gene)
    D['Symbol'].append(M[gene])

    D['Variational'].append(Maps[gene][0])
    D['KNN'].append(Maps[gene][1])
    D['Core'].append(Maps[gene][2])

    D['Cluster ID'].append(l)
    D['Final Core'].append(c)

D = pd.DataFrame(D)
D.to_excel('IBD_12.xlsx')
'''
