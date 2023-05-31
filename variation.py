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


def find_vars(GList, con, dis, V):

    # con = GList[0]
    # dis = GList[1]
    con_max = max([con[u][v]['weight'] for (u, v) in con.edges()])
    dis_max = max([dis[u][v]['weight'] for (u, v) in dis.edges()])

    for (u, v) in con.edges():
        con[u][v]['weight'] = con[u][v]['weight'] / float(con_max)

    for (u, v) in dis.edges():
        dis[u][v]['weight'] = dis[u][v]['weight'] / float(dis_max)

    g_names = list(con.nodes())
    D = {}
    for u in g_names:
        uC = list([con[u][gene]['weight'] for gene in g_names if u != gene]) + \
             list([con[gene][u]['weight'] for gene in g_names if u != gene])

        uD = list([dis[u][gene]['weight'] for gene in g_names if u != gene]) + \
             list([dis[gene][u]['weight'] for gene in g_names if u != gene])

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

