import numpy as np
import networkx as nx
import random
import mygene

from sklearn.metrics import mean_squared_error, mean_absolute_error
from girvan import *
from scipy import sparse
from copy import deepcopy
# from pygenome import sg


def percentile(dicts, per):
    a = list(dicts.values())
    val = np.percentile(a, per)
    # print ('Percentile facts. ', np.mean(a), np.min(a), np.max(a), val)
    return val


def find_vars2(GList, con, dis, V, g_names):
    con = deepcopy(con) / np.max(con)
    dis = deepcopy(dis) / np.max(dis)

    D = {}
    for i in range(len(g_names)):
        gene = g_names[i]

        uC = deepcopy(con[i, :] + con[:, i])
        uD = deepcopy(dis[i, :] + dis[:, i])

        D[gene] = mean_squared_error(uC, uD)

    GList.append(D)

    val = percentile(D, V)
    vars = [key for key in D.keys() if D[key] >= val]
    return vars, GList


def find_vars(GList, con, dis, V):

    con_max = max([con[u][v]['weight'] for (u, v) in con.edges()])
    dis_max = max([dis[u][v]['weight'] for (u, v) in dis.edges()])

    for (u, v) in con.edges():
        con[u][v]['weight'] = con[u][v]['weight'] / float(con_max)

    for (u, v) in dis.edges():
        dis[u][v]['weight'] = dis[u][v]['weight'] / float(dis_max)

    g_names = list(con.nodes())
    D = {}
    for u in g_names:
        # uC = list([con[u][gene]['weight'] for gene in g_names if u != gene and con.has_edge(u, gene)])
        #
        #      list([con[gene][u]['weight'] for gene in g_names if u != gene])
        #
        # uD = list([dis[u][gene]['weight'] for gene in g_names if u != gene]) + \
        #      list([dis[gene][u]['weight'] for gene in g_names if u != gene])

        ind, out = [], []
        for gene in g_names:
            if con.has_edge(u, gene):
                out.append(con[u][gene]['weight'])
            else:
                out.append(0.0)

            if con.has_edge(gene, u):
                ind.append(con[gene][u]['weight'])
            else:
                ind.append(0.0)

        uC = deepcopy(out + ind)

        ind, out = [], []
        for gene in g_names:
            if dis.has_edge(u, gene):
                out.append(dis[u][gene]['weight'])
            else:
                out.append(0.0)

            if dis.has_edge(gene, u):
                ind.append(dis[gene][u]['weight'])
            else:
                ind.append(0.0)

        uD = deepcopy(out + ind)

        D[u] = mean_squared_error(uC, uD)

    GList.append(D)

    val = percentile(D, V)
    vars = [key for key in D.keys() if D[key] >= val]
    return vars, GList


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

    print ([len([labels[gene] for gene in labels.keys() if labels[gene] == i]) for i in range(c)])

    Z = nx.DiGraph()
    for gene in G.nodes():
        Z.add_node(gene, weight = labels[gene])
    Z.add_edges_from(list(G.edges()))

    return labels, c, Z


def convert(gene):
    gene = gene.upper()

    try:
        desc = sg.stdgene[gene]
        desc = desc.short_description
    except:
        desc = ''

    return [gene, desc]


def cores(H, out, ind):
    I = nx.DiGraph()
    I.add_nodes_from(list(H.nodes()))
    I.add_edges_from(list(H.edges()))

    while True:
        l = []
        flag = True
        for u in I.nodes():
            if I.out_degree(u) < out or I.in_degree(u) < ind:
                l.append(u)
                flag = False

        I.remove_nodes_from(l)
        if flag:
            break

    return I


def dcore_main2(G):
    # i, j := in-degree, out-degree
    imin, jmin = min([G.in_degree(u) for u in G.nodes()]), min([G.out_degree(u) for u in G.nodes()])
    imax, jmax = max([G.in_degree(u) for u in G.nodes()]), max([G.out_degree(u) for u in G.nodes()])
    # print (imin, imax, jmin, jmax)

    in_max = [imin, jmin, None]
    out_max = [imin, jmin, None]

    for i in range(imin, imax + 1):
        for j in range(jmin, jmax + 1):
            H = cores(G, j, i)

            if len(H) > 0 and len(H) < len(G):
                if (in_max is None) or (i > in_max[0] and j >= in_max[1]):
                    in_max = [i, j, deepcopy(H)]

                if (out_max is None) or (i >= out_max[0] and j > out_max[1]):
                    out_max = [i, j, deepcopy(H)]

    if (in_max[2] is not None) and ((in_max[0] >= out_max[0] and in_max[1] >= out_max[1]) or out_max[2] is None):
        return in_max

    return out_max

