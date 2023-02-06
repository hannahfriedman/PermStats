from compress import *
from rsk import count_inc_seq
from random_walk import rep_w_ij_kl
from perm_stats import w_ij_kl
import numpy as np
from permutation import Permutation
from itertools import permutations
from math import factorial
from compress import compress


def inc_seq_in_supp(i, j, k, l, n):
    perms = []
    for sigma in Permutation.group(n):
        if sigma(k) == i and sigma(l) == j and count_inc_seq(sigma, 2, n) > 0:
            perms.append(sigma)
    return perms

i = 5
j = 6
k = 1
l = 3
n = 6
M = np.zeros((n*(n - 1), n*(n - 1)))
# for k in range(1, n+1):
#     for l in range(1, n + 1):
#         if k != l:
#             perms1 = inc_seq_in_supp(i, j, k, l, n)
#             if len(perms1) <= 10:
#                 print(k, l, len(perms1))
perms2 = inc_seq_in_supp(i, j, 5, 6, n)
# # print([perm for perm in perms2 if perm not in perms1])
# coeff = [1, 1, -2]
# shift = [perm * Permutation(1, 2, 3, 6, 4, 5) for perm in perms1]
# print(len(perms2), len(shift))
# print([ i for i in perms2 if i not in shift])
# for perm in shift:
#     if count_inc_seq(perm, 1, n) == 0:
#         print(perm, count_inc_seq(perm, 2, n))
for i in range(len(perms2)):
    M += (perms2[i].sign  +1/23) * compress(perms2[i], 2, n)
for m in M:
    print(list(np.round_(m)))
