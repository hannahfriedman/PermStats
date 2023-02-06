import numpy as np
from permutation import Permutation
from dft import inverse_dft
from rsk import count_inc_seq

n = 4
coeff = inverse_dft([
    np.array([[0]]), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]),
    np.array([[1, 0], [0, 0]]), np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]]), np.array([[0]])], 4)
i = 0
for sigma in Permutation.group(n):
    if count_inc_seq(sigma, 2, n) == 0:
        print(24 * coeff[i], sigma)
    i += 1
