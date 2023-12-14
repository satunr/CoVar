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
import community as community_louvain

from networkx.algorithms.community.centrality import girvan_newman
from Nearest import *
from variation import *
from scipy import stats
from constant import *
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


def find_cores(G, cut_off):
    G.remove_edges_from([(u, v) for (u, v) in G.edges() if G[u][v]['weight'] < cut_off])

    # I = list(nx.isolates(G))
    # G.remove_nodes_from(I)

    labels, lC, G = comm2(G, None)
    Cores = [[gene for gene in labels.keys() if labels[gene] == i] for i in range(lC)]

    return G, Cores, labels


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


D_all = pickle.load(open('Meta_Module3.p', 'rb'))
[D0, E0] = deepcopy(D_all[0])
D0 = {str(key[1]): D0[key] for key in D0.keys()}
print (D0)

H = nx.read_gml('Meta3.gml')
print (H.nodes())
partition = community_louvain.best_partition(H)
print (partition)

how_many_runs = 25
# mg = mygene.MyGeneInfo()

# Combined analysis
CombinedKnn = nx.DiGraph()
Maps = {}
C = []
for ind in range(how_many_runs):
    print ('Iteration', ind)

    # < Control, Disease, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >
    # B = pickle.load(open('GList-Macro-a-' + str(ind) + '.p', 'rb'))
    B = pickle.load(open('Run2-' + str(ind) + '.p', 'rb'))
    [_, _, gnames, _, Va, Knn, _, C] = deepcopy(B)
    print (len(Knn.nodes()), len(Knn.edges()))
    # exit(1)

    # M = maps(gnames)

    Cores = []
    for gene_list in C:
        # print (gene_list)
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

    nx.write_gml(CombinedKnn, 'Comb_Meta_App2.gml')
    pickle.dump(Maps, open('Maps_App2.p', 'wb'))


# Remap nodes (Ensemble to symbol)
# --------------------------------
'''
M = maps(gnames)
pickle.dump(M, open('ens_symbol.p', 'wb'))
'''

# Create Combined network
# -----------------------
Maps = pickle.load(open('Maps_App2.p', 'rb'))
gnames = list(Maps.keys())

# M = pickle.load(open('ens_symbol.p', 'rb'))

G = nx.read_gml('Comb_Meta_App2.gml')
# G = nx.relabel_nodes(G, {gene: M[gene] for gene in gnames})
print (list(G.nodes()))

# Find cores in combined network
G, Cores, labels = find_cores(G, confidence)
print (len(G), len(G.edges()))
print ([G.nodes[u]['weight'] for u in G.nodes()])
nx.write_gml(G, 'Comb_Meta_Trimmed_2.gml')
# input('')

condense = []
for each in Cores:
    for gene in each:
        condense.append(gene)

print ([len(c) for c in Cores])

# Create spreadsheet
# ------------------

D = {'Ensemble': [], 'Symbol': [], 'KNN': [], 'Variational': [],
     'Core': [], 'Final Core': [], 'CoVar Cluster ID': [], 'Meta Cluster ID': []}

for gene in Maps.keys():
    # if str(gene) not in H.nodes():
    #     D['Meta Cluster ID'].append(-1)
    # else:
    #     D['Meta Cluster ID'].append(partition[str(gene)])

    if gene in D0.keys():
        D['Meta Cluster ID'].append(D0[gene])
    else:
        D['Meta Cluster ID'].append(-1)

    if gene in labels.keys():
        # l = labels[M[gene]]
        l = labels[gene]
    else:
        l = -1

    if gene in condense:
        c = 1
    else:
        c = 0

    D['Ensemble'].append(int(gene))
    D['Symbol'].append(int(gene))

    D['Variational'].append(Maps[gene][0])
    D['KNN'].append(Maps[gene][1])
    D['Core'].append(Maps[gene][2])

    D['CoVar Cluster ID'].append(l)
    D['Final Core'].append(c)

D = pd.DataFrame(D)
D.to_excel('Meta_3_App2.xlsx')
