import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


def cv_find(l):
    if float(np.mean(l)) == 0:
        return 0.0

    return float(np.std(l)) / float(np.mean(l))


'''
f = open('countMatrix2.txt')

isFirst = True
M = []
for line in f.readlines():
    if isFirst:
        isFirst = False
        continue

    m = [float(each) for each in line.split('\t')]
    M.append(m)

M = np.array(M)
M = deepcopy(M[:6, :]).T
print (np.shape(M))

Y1, Y2 = [], []
for i in range(M.shape[0]):
    Y1.append(cv_find(M[i, :3]))
    Y2.append(cv_find(M[i, 3: 6]))

Y1 = np.histogram(Y1, bins = [float(i) / 10.0 for i in range(11)])[0]
Y2 = np.histogram(Y2, bins = [float(i) / 10.0 for i in range(11)])[0]
print (Y1)
print (Y2)

plt.bar([i * 0.1 - 0.0125 for i in range(10)], [each / float(sum(Y1)) for each in Y1], width = 0.025,
        label = 'Control')
plt.bar([i * 0.1 + 0.0125 for i in range(10)], [each / float(sum(Y2)) for each in Y2], width = 0.025,
        label = 'Perturbed')

plt.xticks([i / 10.0 for i in range(10)], [str((i + 1)/ 10.0) for i in range(10)])
plt.xlabel('Interval of CV')
plt.ylabel('Frequency')
plt.legend()
plt.tight_layout()
plt.savefig('CV.png', dpi = 300)
plt.show()
'''

# Correlation within samples
width = 0.1
C = [(0.990391072020618, 0.0), (0.9916279451614871, 0.0), (0.9996610823043873, 0.0)]
P = [(0.9997097947357286, 0.0), (0.9996100037835742, 0.0), (0.9994828096113242, 0.0)]

plt.bar([i - width / 2.0 for i in range(len(C))], [val[0] for val in C], width, label = 'Control', color = 'C0')
plt.bar([i + width / 2.0 for i in range(len(P))], [val[0] for val in P], width, label = 'Perturbed', color = 'C1')
plt.xticks([0, 1, 2], [(0, 1), (0, 2), (1, 2)])
plt.xlabel('Sample pairs')
plt.ylabel('Pearson correlation coefficient')

plt.ylim([0.95, 1.01])

plt.legend()
plt.tight_layout()
plt.savefig('corr.png', dpi = 300)
plt.show()

