# CoVar: Network Analysis Tool

__authors__ = "Satyaki Roy, Shehzad Sheikh, and Terrence Furey"

__copyright__ = "Furey and Sheikh Labs, UNC Chapel Hill USA"

__email__ = "satyakir@unc.edu, shehzad_sheikh@med.unc.edu, and tsfurey@email.unc.edu"


----------------------------------------------------------------------------------------------------------------------------------------



## Overview.
<p align="justify"> CoVar is a network analysis tool that leverages machine learning and network inference on genomic data, such as gene expression data. It employs the GENIE3 [1] machine learning-based network inference approach to construct a directed network, capturing regulatory interactions between genes. The tool identifies variational genes—genes showing differences in network properties between control and perturbed samples. CoVar then establishes a nearest-neighbor network of variational genes and their strongest interacting neighbors. Within this network, it defines core genes, characterized by both coordination (strong mutual interactions) and reachability (regulatory paths to nearest neighbor network genes). </p>



----------------------------------------------------------------------------------------------------------------------------------------





## Dataset (Input / Output) and Scripts.

Input is a text file with delimited by tabs. The rows are samples and the columns are genes. The rows contain both the control and perturbed samples. 

1. **constant.py:**     It has the list of parameters

2. **Main.py:**         Main script that find the variational and core genes in the different runs. The results are saved in the form of a pickle file of lists in                           the following format *<Control matrix, Disease matrix, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >*. File                                   "GList-i.p" is the result from the i-th run.

3. **read.py:**         It reads and preprocesses the expression dataset.

4. **Genie.py:**        It runs the GENIE3 model to find the inter-gene interaction.

5. **variation.py:**    It has the modules for (1) variational gene identification and (2) network community detection.

6. **Nearest.py:**      It runs the nearest neighbor module.

7. **Combined.py:**     Given the individual sample runs, it finds the nearest neighbor and core genes based on the combination of the runs.


----------------------------------------------------------------------------------------------------------------------------------------





## Modules
There are two modules:

### Identification of variational and core genes (Main.py)

This code eliminates lowly expressed genes, invokes the GENIE3 module to identify, separately on the control and perturbed datasets, the weights w[u, v] representing the influence of each gene u on gene v. It then (1) marks the genes with the highest variation in network characteristics between the control and perturbed networks, (2) finds a nearest neighbor network comprising the variational genes and their nearest neighbors, and (3) identifies strongly connected clusters (or communities) within the neartest neighbor networks and core genes within each community.

The above steps are followed in each run. As mentioned earlier, the results are saved in the form of a pickle file of lists in the following format *<Control matrix, Disease matrix, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >*. File "GList-i.p" is the result from the i-th run.

### Combined analysis (Combined.py)

To make a robust inference, this module finds an aggregate nearest neighbor network across the 25 runs. It then finds communities within this network and core genes within each community. It finally produces a spreadsheet with the following information on the nearest neighbor genes: 

['Gene symbol', 'Cluster ID', 'Gene description', 'Frequency of that gene as core in 25 runs', 'Frequency of that gene as nearest neighbor gene in 25 runs', 
 'In-degree', 'Out-degree', 'Frequency of that gene as variational in 25 runs', 'Final core Gene']

----------------------------------------------------------------------------------------------------------------------------------------

At present the code employs the yeast dataset (countMatrix2.txt) from the study on blocking the oxidation-phosphorylation pathways [2]. 

----------------------------------------------------------------------------------------------------------------------------------------


[1] V. Huynh-Thu, et al. "Inferring regulatory networks from expression data using tree-based methods." PloS one 5.9 (2010): e12776.

[2] S. Liu S, et al. OXPHOS deficiency activates global adaptation pathways to maintain mitochondrial membrane potential. EMBO Rep 2021 Apr 7;22(4):e51606. PMID: 33655635
