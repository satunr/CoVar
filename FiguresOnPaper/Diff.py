import numpy as np
import pandas as pd
import pickle
import networkx as nx
import matplotlib.pyplot as plt

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from copy import deepcopy
from scipy.stats import pearsonr
from scipy import stats


def readf(fname):
    f = open(fname, 'r')
    l = f.readlines()

    g_names = l[0]
    g_names = [gene.replace('\n', '') for gene in g_names.split()]

    A = []
    for i in range(1, len(l)):
        A.append([float(val) for val in l[i].split()])

    A = np.array(A)

    return A, g_names


def pearson_find(X, Y):
    p, _ = pearsonr(X, Y)
    return float(p)


def nearest_neighbor(coor_control, coor_disease):
    common = {}

    for u in range(coor_control.shape[0]):

        dc = {v: pearson_find(coor_control[u, :], coor_control[v, :]) for v in range(coor_control.shape[0])}
        dd = {v: pearson_find(coor_disease[u, :], coor_disease[v, :]) for v in range(coor_control.shape[0])}
        common[str(u)] = np.mean(list({v: abs(dd[v] - dc[v]) for v in range(coor_control.shape[0])}.values()))

    return common


def write(Prob_randomness, Prob_randomness_0, fname):

    D = {'Gene Set': [], 'Control to random score': [],
         'Perturbed to random score': [], 'Absolute difference of scores': []}

    for key in Prob_randomness.keys():
        sc = abs(Prob_randomness_0[key][0] - Prob_randomness[key][0])
        gs = [str(each) for each in Prob_randomness_0[key][1]]

        D['Gene Set'].append(', '.join(gs))
        D['Control to random score'].append(Prob_randomness_0[key][0])
        D['Perturbed to random score'].append(Prob_randomness[key][0])
        D['Absolute difference of scores'].append(sc)

    D = pd.DataFrame(D)
    D.to_excel(fname)


def hierarchical_clustering_best_config(data, max_clusters = 199):
    """
    Perform hierarchical clustering and choose the best configuration based on silhouette score.

    Parameters:
    - gene_expression_data: Gene expression data as a numpy array.
    - max_clusters: Maximum number of clusters to consider.

    Returns:
    - Best clustering configuration (number of clusters) and the corresponding labels.
    """
    best_score = -1
    best_config = 0
    best_labels = None

    for n_clusters in range(2, max_clusters + 1):
        # Perform hierarchical clustering
        model = AgglomerativeClustering(n_clusters=n_clusters)
        labels = model.fit_predict(data)

        # Calculate silhouette score
        score = silhouette_score(data, labels)
        # print (n_clusters, score)

        # Update best configuration if the current score is higher
        if score > best_score:
            best_score = score
            best_config = n_clusters
            best_labels = labels

    return best_config, best_labels


def find_t(gene_expression_data, indices):

    PCC = []
    for i in range(len(indices) - 1):
        for j in range(i + 1, len(indices)):
            d1 = gene_expression_data[indices[i]]
            d2 = gene_expression_data[indices[j]]
            corr, _ = pearsonr(d1, d2)

            PCC.append(corr)

    t = np.mean(PCC) / (np.std(PCC) / np.sqrt(len(PCC)))
    return t


def each_dataset(data, best_config, best_labels):

    Prob_randomness = {}
    for c in sorted(range(best_config)):
        indices_cluster = [i for i in range(len(best_labels)) if best_labels[i] == c]
        if len(indices_cluster) < 3:
            continue

        t_cluster = find_t(data, indices_cluster)

        prob_randomness = 0
        for j in range(how_many_random):
            indices_random = np.random.choice([i for i in range(n_genes)], size = len(indices_cluster)).tolist()
            t_random = find_t(data, indices_random)

            if t_random > t_cluster:
                prob_randomness += 1.0 / float(how_many_random)

        # print (len(indices_cluster), t_cluster, prob_randomness)
        Prob_randomness[c] = [prob_randomness, indices_cluster]

    return Prob_randomness


# This code generates a scatter plot showing the 
# differential coexpression of a gene against its 
# CoVar variationality (measured in terms of mean 
# squared error in weights in the GENIE networks generated 
# from the control and perturbed datasets)

fname = 'Meta3.txt'
gene_expression_data, _ = readf(fname)
gene_expression_data = deepcopy(gene_expression_data.T)
print (gene_expression_data.shape)

control = deepcopy(gene_expression_data[:, :15])
perturb = deepcopy(gene_expression_data[:, 15:])

# Mean correlation
# ================
common = nearest_neighbor(control, perturb)
print (common)

X, Y = [], []
for r in range(1):
    # D = pickle.load(open('/Users/satyakiroy/PycharmProjects/pythonProject/Genomics/CoVar2_Run/Run2-' + str(r) + '.p', 'rb'))
    D = pickle.load(open('/Users/satyakiroy/PycharmProjects/pythonProject/Genomics/CoVar2_Run/IBD-0.p', 'rb'))
    MSE = D[3]
    print (MSE)

    x = [MSE[gene] for gene in sorted(MSE.keys())]
    y = [common[gene] for gene in sorted(MSE.keys())]
    plt.scatter(x, y, marker = 'x', s = 1)
    X.extend(x)
    Y.extend(y)

plt.xlabel('Mean square error')
plt.ylabel('Differential correlation')

pcc = round(pearson_find(X, Y), 3)
print (pcc)

plt.title('Pearson correlation coefficient ' + str(pcc))
plt.tight_layout()
plt.show()
