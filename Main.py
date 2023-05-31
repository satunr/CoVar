import networkx as nx
import numpy as np
import random
import matplotlib.pyplot as plt
import collections
import stat
import pickle
import math
import pandas as pd

from networkx.algorithms.community.centrality import girvan_newman
from Nearest import *
from variation import *
from scipy import stats
from copy import deepcopy
from itertools import permutations, combinations
from Genie import GENIE3, create_graph
from constant import *
from read import *

np.printoptions(suppress = True, precision = 3)


def color_nodes(Knn, Cr, labels):
    print ('***', list(set(labels.values())))

    H = nx.DiGraph()
    D = {'Gene': [], 'Out-degree': [], 'In-degree': [], 'Core': [], 'Labels': []}

    for each in Cr:
        print (each)

        for u in each:
            H.add_node(u, weight = 1)

            D['Gene'].append(u)
            D['In-degree'].append(Knn.in_degree(u))
            D['Out-degree'].append(Knn.out_degree(u))
            D['Core'].append(1)
            D['Labels'].append(labels[u])

    L = [u for u in list(Knn.nodes()) if u not in H.nodes()]
    for u in L:
        H.add_node(u, weight = 0)
        D['Gene'].append(u)
        D['In-degree'].append(Knn.in_degree(u))
        D['Out-degree'].append(Knn.out_degree(u))
        D['Core'].append(0)
        D['Labels'].append(labels[u])

    H.add_edges_from(list(Knn.edges()))
    D = pd.DataFrame(D)
    # D.to_excel('Res_IBD.xlsx')
    return H


def main(ind):
    # STEP 1: Read count matrix
    data, g_names = readf(fname)

    # Extract necessary columns (or samples)
    extract = [i for i in range(indices[0][0], indices[0][1])] + [i for i in range(indices[1][0], indices[1][1])]
    print (extract)

    data = data[extract, :]
    print (np.shape(data))

    data, g_names = remove_noncoding(data, g_names)
    print (len(g_names))

    # Find TPM
    # lengths = read_annotations('annotationMatrix2.csv')
    # data = find_tpm(lengths, data, g_names)
    # exit(1)

    # Find gene symbols
    m, g_names = mappings(g_names)
    pickle.dump([g_names, m], open('m.p', 'wb'))

    [g_names, m] = pickle.load(open('m.p', 'rb'))
    g_names = [m[gene] for gene in g_names]
    # check(g_names)

    # STEP 2: Eliminate lowly expressed genes
    D = None
    # D = pickle.load(open('M-values.p', 'rb'))
    data, g_names, D = eliminate(data, g_names, cut_off, cv_cutoff, D)
    pickle.dump([data, g_names], open('Y.p', 'wb'))
    print (len(g_names))

    [data, g_names] = pickle.load(open('Y.p', 'rb'))
    print ('Number of genes left. ', len(g_names))
    print (len(g_names), np.shape(data))

    # STEP 3: Normalize dataset
    data[-1, :] += 0.01
    noise = np.random.normal(0, 0.0000001, data.shape)
    data = data + noise
    data = data/np.max(data)

    # STEP 3: Store control and disease data (gene names, normalized count matrix, and genie network)
    S = [[], [], []]

    # < Control, Disease, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >
    GList = []
    con, dis = None, None
    for i in range(2):
        S[0].append(g_names)
        S[1].append(data[X[i]: X[i + 1], :])

        # STEP 4a: Create genie matrix and isolate gene names
        VIM3 = GENIE3(data[X[i]: X[i + 1], :], nthreads = nth, tree_method=tree_method, ntrees = ntrees)
        VIM3 = VIM3/np.max(VIM3)
        print (np.shape(VIM3))

        GList.append(VIM3)

        # STEP 4b: Create genie graph
        G = create_graph(VIM3, g_names)
        S[2].append([G, VIM3])

        if i == 0:
            con = deepcopy(G)
        elif i == 1:
            dis = deepcopy(G)
            GList.append(G)

    vg, GList = find_vars(GList, con, dis, V)
    GList.append(vg)

    # STEP 5: Find k-nearest neighbor graph on disease graph
    knn = nearest2(dis, vg, K, approach)
    print (len(knn), len(knn.edges()))
    GList.append(knn)

    # Find communities
    labels, lC, _ = comm2(knn, GList)
    GList.append(labels)

    # STEP 6: Find cores in each community
    Cr = []
    for i in range(lC):
        subset = [j for j in knn.nodes() if labels[j] == i]
        if len(subset) < 2:
            c = deepcopy([subset[0][:]])
        else:
            c = deepcopy(list(nx.k_core(create(knn, subset)).nodes()))
        Cr.append(c)
    print ([len(Cr[i]) for i in range(len(Cr))])

    GList.append(Cr)
    pickle.dump(GList, open('GList-Macro-a-' + str(ind) + '.p', 'wb'))


for ind in range(1):
    # Each run
    print ('Run.', ind)
    main(ind)

    '''
    [_, m] = pickle.load(open('m.p', 'rb'))
    m = {gene[4:]: m[gene] for gene in m.keys()}
    print (m)
    '''

    '''
    # GList <Control, Perturbed, g_names, MSE, variational, KNN, Labels, Core>
    GList = pickle.load(open('GList-Macro-a-' + str(ind) + '.p', 'rb'))
    Knn = deepcopy(GList[5])
    labels = deepcopy(GList[6])

    print (len(Knn), len(Knn.edges()))
    print (list(Knn.nodes()))

    Cr = deepcopy(GList[7])

    Knn = color_nodes(Knn, Cr, labels)
    # nx.write_gml(Knn, 'Knn.gml')
    '''

