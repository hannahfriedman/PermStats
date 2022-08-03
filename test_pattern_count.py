from perm_stats import count_occurrences
from perm_stats import length
from permutation import Permutation
from dft import dft
import numpy as np

# n = 5
# for sigma in Permutation.group(n):
#     print(count_occurrences(Permutation(3, 2, 1), n, 3)(sigma, n), length(sigma, n))

n = 6
for tau in Permutation.group(3):
    print(tau)
    print(np.round_(np.linalg.eig(dft(count_occurrences(tau, n, 3), n)[7])[0], 10))
