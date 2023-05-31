import numpy as np
import pickle, math
import pandas as pd
import mygene
import time
import ast

# from pygenome import sg
from scipy.stats import *


def update_(X, gene, mv, s, val):
    X['Gene'].append(gene)
    X['M-value'].append(mv)
    X['Sum(read count)'].append(s)
    X['Include'].append(val)
    return X


def eliminate(A, g_names, cut_off, cv_cutoff, D):

    X = {'Gene': [], 'M-value': [],
         'Sum(read count)': [], 'Include': []}

    if D is None:
        # Calculate M-value
        D = mvalue(A, g_names)
        print(D)
        pickle.dump(D, open('M-values.p', 'wb'))

    # Remove genes and samples with all 0 expression
    rem_gene, elim = [], []
    for i in range(len(g_names)):
        mv = D.iloc[i]['M-value']

        if sum(A[:, i]) > cut_off and mv > cv_cutoff:
            rem_gene.append(i)
            X = update_(X, g_names[i], mv, sum(A[:, i]), 1)
        else:
            elim.append(i)
            X = update_(X, g_names[i], mv, sum(A[:, i]), 0)

    A = A[:, rem_gene]
    g_names = remove_genes(g_names, rem_gene)

    X = pd.DataFrame(X)
    X.to_excel('Filter.xlsx')

    return A, g_names, D


def readf(fname):
    f = open(fname, 'r')
    l = f.readlines()

    g_names = l[0]
    g_names = [gene.replace('\n', '') for gene in g_names.split()]
    # print (g_names)

    # X = [('MATR3', 'ENSG00000280987.4'), ('PRORP', 'ENSG00000258790.1'), ('ST6GALNAC6', 'ENSG00000257524.6')]
    # X = [g_names.index(x[1]) for x in X]

    A = []
    for i in range(1, len(l)):
        A.append([float(val) for val in l[i].split()])

    A = np.array(A)
    # print (len(g_names), A.shape)
    # print (X)

    return A, g_names


def mvalue(data, g_names, c = 0.0000001):
    nSamples = np.shape(data)[0]
    D = {'Gene': [], 'M-value': []}

    for i in range(len(g_names)):
        print (float(i) / float(len(g_names)), g_names[i])

        gene = g_names[i]
        m = np.mean(data[:, i])

        if m > 0:
            D['Gene'].append(gene)

            V = []
            for k in range(len(g_names)):
                if i == k or np.mean(data[:, k]) == 0:
                    continue

                v = [math.log((float(data[j, i]) + c) / (float(data[j, k]) + c), 2) for j in range(nSamples)]
                v = np.std(v)
                V.append(v)
            D['M-value'].append(float(sum(V)) / float(len(V)))

        else:
            D['Gene'].append(gene)
            D['M-value'].append(0)

    return pd.DataFrame(D)


def remove_genes(g_names, rem):
    g1_names = []
    for i in range(len(g_names)):
        if i in rem:
            g1_names.append(g_names[i])
    return g1_names


def mappings(glist):
    ginfos = {}
    mg = mygene.MyGeneInfo()

    # g_names = [gene[:gene.index('.')] for gene in glist]
    g_names = [gene[5:] for gene in glist]
    # ginfo = mg.querymany(g_names, scopes='ensembl.gene', returnall=True, species='human')

    ginfo = mg.querymany(g_names, scopes='ensembl.gene', returnall=True)
    ginfo = ginfo['out']
    how_many = [0, 0]

    for i in range(len(ginfo)):

        if 'symbol' in ginfo[i].keys():
            ginfos[ginfo[i]['query']] = ginfo[i]['symbol']
            how_many = [how_many[0] + 1, how_many[1]]
        else:
            ginfos[ginfo[i]['query']] = ginfo[i]['query']

    return ginfos, g_names


def read_annotations(fname):
    D = pd.read_csv(fname)
    lengths = {}
    for i in range(len(D)):
        lengths[D.loc[i]['GeneID']] = D.loc[i]['Length']
    return lengths


def remove_noncoding(data, g_names):
    pc_genes = [i for i in range(len(g_names)) if 'gene-Y' in g_names[i]]
    return data[:, pc_genes], [g_names[i] for i in pc_genes]


def find_tpm(lengths, data, g_names, factor = 1000000):
    data = data.astype(float)

    # Step 1: Scale by gene length
    for i in range(np.shape(data)[1]):
        data[:, i] = data[:, i] / float(lengths[g_names[i]])
    # print ('Step 1:\n', data)

    # Step 2: Scale by sequencing depth
    depth = {i: np.sum(data[i, :]) for i in range(np.shape(data)[0])}
    # print (depth)

    for i in range(np.shape(data)[0]):
        data[i, :] = data[i, :] / (float(depth[i]) / factor)
    # print (data)
    # print ('Step 2:\n', [np.sum(data[i, :]) for i in range(np.shape(data)[0])])
    return data


def check(g_names):
    cnt = {}
    for gene in g_names:
        if gene not in cnt.keys():
            cnt[gene] = 0
        cnt[gene] = cnt[gene] + 1
    cnt = {gene: cnt[gene] for gene in cnt.keys() if cnt[gene] > 1}
    # print (len(cnt), sorted(cnt.items(), key=lambda x: x[1], reverse=True))
