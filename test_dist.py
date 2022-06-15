from dft import dft
from dft import dft_natural
from dft import sort_partitions
from dft import generate_partitions
from dft import dft_matrix
from dft import one_local_dft_natural
from dft import projection
from dft import inverse_dft
from dft import inverse_dft_natural
from permutation import Permutation
from random_walk import representation
from perm_stats import length
from perm_stats import distance_from_standard
from perm_stats import w_ij_kl
from random_walk import w_ij_mat
from random_walk import w_ij_kl_mat
from tableau import Tableau
import misc
import numpy  as np
from math import factorial
import matplotlib.pyplot as plt
from misc import function_to_vector
from misc import distribution
from itertools import combinations
from sympy import Matrix
from decompress import T_w_ij_kl
n = 5
sn = [g for g in Permutation.group(n)]
def compare_with_action(set1, set2, sigma, tau):
    set3 = set([sigma * pi * tau for pi in set2])
    return set(set1) == set3
def counts(vector, val):
    return [sn[i] for i in range(len(sn)) if vector[i] == val]

# print(projection(distance_from_standard([3, 1, 1], n), n, 1))
# Parititons for n = 4:
# [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]


# for partition in sort_partitions(generate_partitions(n)):
# for partition in [(3, 1, 1), (2, 2, 1)]:
ft1 = dft_natural(distance_from_standard((3, 1, 1), n), n)
ft2 = dft_natural(distance_from_standard((2, 2, 1), n), n)
for m in ft1:
    print(m)
for m in dft_natural(misc.vector_to_function(inverse_dft_natural(ft1, n), n), n):
    print(m)

diff = [misc.similar(ft1[i], ft2[i]) for i in range(len(ft1))]
preimage = inverse_dft_natural(diff, n)

int_pre_image = np.round_(360*np.array(preimage), 5)
two_local_part = projection(misc.vector_to_function(int_pre_image, n), n, 0, 1, 2, 3)
print(two_local_part)
neg_sum = 0
pos_sum = 0
for i in range(len(two_local_part)):
    if two_local_part[i] <= 0:
        neg_sum += two_local_part[i]
    else:
        pos_sum += two_local_part[i]
print(neg_sum, pos_sum)

# dft_vec = T_w_ij_kl(n)@ np.array(two_local_part)
# print(dft_vec)
# total = np.array(np.zeros((1, factorial(n))))
# count = 0
# for i in range(1, n+1):
#     for j in range(1, n+1):
#         if i != j:
#             for k in range(1, n+1):
#                 for l in range(1, n+1):
#                     if k != l:
#                         total += dft_vec[count] * np.array(function_to_vector(w_ij_kl(i, j, k, l), n))
#                         count += 1
# print(total)
# for n in range(4, 5):
#     m = Matrix(T_w_ij_kl(n))
#     print(m)
#     new, blurg = m.rref()
#     m = np.array(new).astype(np.float64)
#     print(np.linalg.matrix_rank(m))
#     print(m[:23])

# diff2 = [misc.similar(ft2[i], ft1[i]) for i in range(len(ft1))]
# preimage2 = inverse_dft_natural(diff2, n)
# print(np.round_(360*np.array(preimage2), 5))
# image = [np.zeros((n-1, n-1)), np.zeros((5, 5)), P1 @ np.linalg.inv(P2),  np.zeros((5, 5)), np.zeros((n-1, n-1)), np.zeros((1, 1))]
# print(inverse_dft(image, n))
    # print(np.round_(np.linalg.eig(dft_natural(distance_from_standard(partition, n), n)[3])[0], 5))
    # print(np.round_(np.linalg.eig(dft_natural(distance_from_standard(partition, n), n)[3])[1], 5))
        # print(np.round_(np.linalg.eig(m)[0], 5))
    # plt.matshow(representation(length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False) - representation(distance_from_standard(partition, n), w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False), cmap='Blues')
    # plt.matshow(representation(length, w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False), cmap='Blues')
    # print(representation(distance_from_standard(partition, n), w_ij_kl_mat, misc.falling_factorial(n, n-2), n, False))
    # print(partition, '\n', representation(distance_from_standard(partition, n), w_ij_mat, n, n, False), '\n',np.round_(np.linalg.eig( representation(distance_from_standard(partition, n), w_ij_mat, n, n, False))[0], 3))
#    print(partition, '\n', np.round_(np.linalg.eig( representation(distance_from_standard(partition, n), w_ij_kl_mat, n*(n-1), n, False))[0], 3))    

# for a, b in combinations(sort_partitions(generate_partitions(n)), 2):
#     print(a, b, distribution(function_to_vector(distance_from_standard(a, n), n)) == distribution(function_to_vector(distance_from_standard(b, n), n)))

# two_two_one = np.array(function_to_vector(distance_from_standard([2, 2, 1], n), n))
# three_one_one = np.array(function_to_vector(distance_from_standard([3, 1, 1], n), n))
# print(distribution(two_two_one))
# print(distribution(three_one_one))

# for sigma, tau in combinations(sn, 2):
#     if compare_with_action(counts(two_two_one, 0), counts(three_one_one, 0), sigma, tau):
#         print(sigma, tau)

# for i in range(0, 6):
#     set1 = counts(two_two_one, i)
#     set2 = counts(three_one_one, i)
#     print(compare_with_action(set1, set2, Permutation(), Permutation(1, 4, 2, 5, 3)))


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

