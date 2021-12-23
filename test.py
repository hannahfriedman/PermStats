## Find eigenvalues of block representations ##

from dft import dft 
from perm_stats import w_ij
from perm_stats import w_ij_kl
import numpy as np
from scipy.linalg import block_diag
from numpy.linalg import eig

n = 4
r = range(1, 5)
for i in r:
    for j in r:
        if i != j:
            for k in r:
                for l in r:
                    if k != l:
                        A = block_diag(*dft(w_ij_kl(i,j,k,l), n))
                        print(i, j, k, l, np.round_(eig(A)[0], decimals = 5))
# n = 4
# r = range(1, n+1)
# for i in r:
#     for j in r:
#         A = block_diag(*dft(w_ij(i,j), n))
#         print(i, j, np.round_(eig(A)[0], decimals = 3))
# A = np.array([[0, 1, 1, 0, 0, 1], [1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 0, 0], [1, 0, 1, 1, 0, 0]])
# print(A)

# n = 6
# for i in range(1, n+1):
#     transform = (dft(w_ij(i, n), n))
#     d = block_diag(transform[0], transform[1])
#     print(d, d*A)
