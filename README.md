# CoVar: Network Analysis Tool

__authors__ = "Satyaki Roy, Shehzad Sheikh, and Terrence Furey"

__copyright__ = "Furey and Sheikh Labs, UNC Chapel Hill USA"

__email__ = "satyakir@unc.edu, shehzad_sheikh@med.unc.edu, and tsfurey@email.unc.edu"


----------------------------------------------------------------------------------------------------------------------------------------

## 1. Overview
<p align="justify"> CoVar is a network analysis tool that leverages machine learning and network inference on genomic data, such as gene expression data. It employs the GENIE3 [1] machine learning-based network inference approach to construct a directed network, capturing regulatory interactions between genes. The tool identifies variational genes—genes showing differences in network properties between control and perturbed samples. CoVar then establishes a nearest-neighbor network of variational genes and their strongest interacting neighbors. Within this network, it defines core genes, characterized by both coordination (strong mutual interactions) and reachability (regulatory paths to nearest neighbor network genes). </p>

----------------------------------------------------------------------------------------------------------------------------------------
 
## 2. Implementation Steps
<p align="justify"> Running CoVar consists of running two main scripts, Main3.py and Combined.py. The first identifies the variational and core genes for a single integration. The second combines information across multiple iterations to generate consensus variational, nearest, and core genes within consensus modules. The following are specific instructions for executing each step. </p>

### Step (a). Identification of variational and core genes 

<p align="justify"> This script eliminates lowly expressed genes and leverages the GENIE3 module to identify the influence of each gene on others. It then marks genes with the highest variation in network characteristics, forms nearest-neighbor networks, and identifies strongly connected clusters and core genes within each community.  </p>

**Execution Narrative**

**Input.** Ensure the expression data is assigned to the variable 'fname' in constant.py.

**Run Main3.py.** Execute the script to perform the identification process.

**Output.** Save results in "Run`<Approach 1 or 2>`-i.p" files, where i represents the run number.

### Step (b): Generate aggregated CoVar network

<p align="justify"> This script ensures robust inference by finding an aggregate nearest neighbor network across multiple runs. It identifies communities within this network and core genes within each community. The final output is a spreadsheet with detailed information on nearest neighbor genes.</p>

**Execution Narrative**

**Input.** Provide the "Run`<Approach 1 or 2>`-i.p" files from the Main3.py runs.

**Run Combined.py.** Execute the script to integrate information and generate the CoVar network.

**Output.** Obtain `<Comb_Trimmed_2.gml>` (CoVar network) and a spreadsheet enumerating variational, nearest neighbor, and core genes.

<p align="justify">By following these steps, you can effectively run CoVar, extracting valuable insights into variational and core genes across multiple iterations. The combination of Main3.py and Combined.py ensures a comprehensive analysis, providing a deeper understanding of biological networks.</p>


----------------------------------------------------------------------------------------------------------------------------------------

## 3. CoVar Project Structure
<p align="justify"> 
 
### Input Data
- Input data is expected in a tab-delimited text file format.
- Rows represent samples, and columns represent genes.
- The file includes both control and perturbed samples.

### Code Files

#### `constant.py`
- Holds a list of parameters used across the project.

#### `Main3.py`
- Main script responsible for identifying variational and core genes in different runs.
- Results are saved as pickle files containing lists with the format: `<Control matrix, Disease matrix, Gene Names, Mean squared, Variational, Knn, Community Labels, Core >`.
- Result (pickle) files are named "Run`<Approach 1 or 2>`-i.p" from the $i$-th run.

#### `read.py`
- Reads and preprocesses the expression dataset.

#### `Genie.py`
- Executes the GENIE3 model to identify inter-gene interactions.

#### `variation.py`
- Contains modules for (1) variational gene identification and (2) network community detection.

#### `Nearest.py`
- Executes the nearest neighbor module.

#### `Combined.py`
- Given individual sample runs, determine nearest neighbor and core genes based on the combination of runs.</p>


----------------------------------------------------------------------------------------------------------------------------------------
At present the code employs the yeast dataset (**countMatrix2.txt**) from the study on blocking the oxidation-phosphorylation pathways [2] or a simulated expression dataset (**Main3.txt**) [3]. 

----------------------------------------------------------------------------------------------------------------------------------------

## References.

[1] V. Huynh-Thu, et al. "Inferring regulatory networks from expression data using tree-based methods." PloS one 5.9 (2010): e12776.

[2] S. Liu S, et al. OXPHOS deficiency activates global adaptation pathways to maintain mitochondrial membrane potential. EMBO Rep 2021 Apr 7;22(4):e51606. PMID: 33655635.

[3] K. Wang, et al. "Meta-analysis of inter-species liver co-expression networks elucidates traits associated with common human diseases." PLoS computational biology 5.12 (2009): e1000616.

