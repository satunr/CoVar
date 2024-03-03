import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from copy import deepcopy
from scipy.stats import pearsonr
from scipy import stats


def cluster(name, exp, t = 5):
    current = deepcopy(exp)
    similarity_matrix = sim(current)

    # Make the similarity matrix symmetric [NOT NEEDED HERE]
    symmetric_matrix = np.minimum(similarity_matrix, similarity_matrix.T)

    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(sch.distance.squareform(1 - symmetric_matrix), method='ward')

    # Extract the cluster assignments
    labels = sch.fcluster(linkage_matrix, t, criterion='maxclust')  # t is the number of clusters

    # Visualization of the dendrogram and clustered data
    dendrogram = sch.dendrogram(linkage_matrix)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xticks(fontsize = 2, rotation = 90)
    
    plt.xlabel('Genes')
    plt.ylabel('Distance')
    
    # Adjust the xticks fontsize
    # plt.tight_layout()
    plt.savefig(name + '_Cluster.png', dpi = 300)
    plt.show()
    
    # Print the cluster assignments
    # print("Cluster Assignments:", labels)

    return labels


def sim(current):

    n = current.shape[0]
    similarity_matrix = np.ones((n, n))
    for i in range(n - 1):
        for j in range(i + 1, n):
            corr, pvalue = pearsonr(current[i], current[j])
            similarity_matrix[i, j] = similarity_matrix[j, i] = round(corr, 2)
            # print (similarity_matrix[i, j], similarity_matrix[j, i])

    similarity_matrix[np.isnan(similarity_matrix)] = 0
    return similarity_matrix


def readf(fname):
    # Read the text file and convert to a vector representation
    f = open(fname, 'r')
    l = f.readlines()

    g_names = l[0]
    g_names = [gene.replace('\n', '') for gene in g_names.split()]
    # print (g_names)

    A = []
    for i in range(1, len(l)):
        # print (len([float(val) for val in l[i].split()]))
        A.append([float(val) for val in l[i].split()])

    A = np.array(A)
    # print (A.shape)

    return A, g_names


def read(M):
    np.fill_diagonal(M, 1)
    return M


# Input expression data (fname) and output 
# dendrogram showing the clusters based 
# on hierarchical clustering

fname = 'Meta.txt'
gene_expression_data, g_names = readf(fname)
gene_expression_data = deepcopy(gene_expression_data.T)

n_genes = gene_expression_data.shape[0]
m_samples = gene_expression_data.shape[1]

# Do not uncomment the Control and Perturbed lines below at the same time.
# Control
# control = deepcopy(gene_expression_data[:, :15])
# print (control.shape)
# labels_control = cluster('Control', control)
# print (labels_control)

# Perturbed
perturb = deepcopy(gene_expression_data[:, 15:])
print (perturb.shape)
labels_perturb = cluster('Perturb',perturb)
print (labels_perturb)
