# List of system parameters
# -------------------------

# Use Extra-Trees method
tree_method = 'ET'
# tree_method = 'RF'

# Number of trees per ensemble
ntrees = 4000

# Count data
fname = 'countMatrix2.txt'
# fname = 'IBD.txt'

# Indices of control and disease networks
# indices = [[0, 5], [5, 10]]
indices = [[0, 3], [3, 6]]

# Arrange sample IDs
X = [0, len([i for i in range(indices[0][0], indices[0][1])]),
     len([i for i in range(indices[0][0], indices[0][1])] + [i for i in range(indices[1][0], indices[1][1])])]
print (X)

# Nearest neighbor approach
approach = 2

# Percentile of the total number of genes considered variational
if approach == 1:
     V = 97.5
else:
     V = 97.5

# Core mode
modes = 'directed'
level = -1

# Number of threads in genie
nth = 10

index = modes[0].upper() + '-' + str(approach)

# Threshold for minimum expression sum
# cut_off = 7.5 * 10.0
cut_off = 60.0 * 6

# cv_cutoff = 0.1
cv_cutoff = 0.32

code_modes = 2
trim_percentile = 0.0

how_many_runs = 25

# Combined analysis
confidence = 3

# Percentile for KNN
K = 99.95
