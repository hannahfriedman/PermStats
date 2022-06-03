from dft import dft
from permutation import Permutation
from random_walk import representation
from perm_stats import length
from perm_stats import distance_from_standard
from random_walk import w_ij_mat
from random_walk import w_ij_kl_mat
from tableau import Tableau
import misc
import numpy  as np
from math import factorial
import matplotlib.pyplot as plt

# Parititons for n = 4:
# [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]

n = 5
for partition in [(5,), (4, 1), (3, 2), (3, 1, 1), (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]:
    # plt.matshow(representation(length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False) - representation(distance_from_standard(partition, n), w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False), cmap='Blues')
    # plt.matshow(representation(length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False), cmap='Blues')
    # print(representation(distance_from_standard(partition, n), w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False))
    #print(partition, '\n', representation(distance_from_standard(partition, n), w_ij_mat, n, n, False), '\n',np.round_(np.linalg.eig( representation(distance_from_standard(partition, n), w_ij_mat, n, n, False))[0], 3))
    print(partition, '\n', np.round_(np.linalg.eig( representation(distance_from_standard(partition, n), w_ij_kl_mat, n*(n-1), n, False))[0], 3))    
#    plt.show()

# for sigma in Permutation.group(n):
#     print(distance_from_standard(sigma, n), length(sigma, n))

# for n in range(3, 8):
#     base = dft(length, n-1)[0][0,0]
#     print(n, base)
#     m =base*np.ones((n, n))
#     additions = representation(length, w_ij_mat, n, n, False) - m
#     for col in range(additions.shape[1]):
#         if col == 1:
#             print(factorial(n - col - 1), factorial(col))
#         additions[:, col] = additions[:, col]/(factorial(n - col - 1) * factorial(col))
#     print(additions)

