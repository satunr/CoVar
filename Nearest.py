import networkx as nx
import operator
import pickle
import pandas as pd
import collections
# import matplotlib.pyplot as plt

from copy import deepcopy
from variation import *


def create(Knn, lst):
    return deepcopy(Knn.subgraph(lst))


def nearest2(G, I, K, approach, kPer = 100):
    knn = nx.DiGraph()
    knn.add_nodes_from(I)
    E2 = {}

    for (u, v) in G.edges():
        E2[(u, v)] = G[u][v]['weight']

    if approach == 1:
        E_sorted = {e: G[e[0]][e[1]]['weight'] for e in list(G.edges())}
        E_sorted = sorted(E_sorted.items(), key=operator.itemgetter(1), reverse = True)
        E_sorted = [val[0] for val in E_sorted]

        node_index = {u: {} for u in I}
        for e in E_sorted:
            if e[0] in I:
                node_index[e[0]][e] = G[e[0]][e[1]]['weight']
            if e[1] in I:
                node_index[e[1]][e] = G[e[0]][e[1]]['weight']

        for u in I:
            val = percentile(node_index[u], K)
            new_edges = [key for key in node_index[u].keys() if node_index[u][key] >= val]

            # for k in range(min(kPer, len(new_edges))):
            #     knn.add_edge(new_edges[k][0], new_edges[k][1])

            for k in range(len(new_edges)):
                knn.add_edge(new_edges[k][0], new_edges[k][1])

    elif approach == 2:

        val = percentile(E2, K)
        E2 = {key: E2[key] for key in E2.keys() if E2[key] >= val}

        k = 0
        while True:

            sub_edge = {key: E2[key] for key in E2.keys() if key[0] in list(knn.nodes()) or key[1] in list(knn.nodes())}
            if len(sub_edge.keys()) == 0:
                break

            e = max(sub_edge, key = sub_edge.get)
            knn.add_edge(e[0], e[1])
            E2.pop(e, None)
            k = k + 1

    isolated = list(nx.isolates(knn))
    knn.remove_nodes_from(isolated)

    return knn
