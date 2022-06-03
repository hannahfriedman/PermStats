from dft import dft
from dft import dft_natural
from perm_stats import w_ij_kl
from perm_stats import w_ij
from perm_stats import cycle_indicator
from random_walk import representation
from random_walk import w_ij_mat
from perm_stats import excedances
from perm_stats import length
from decompress import T_w_ij
import numpy as np
from misc import nullspace
from math import factorial
# i = 1
# j = 2
# k = 1
# l = 2
# n = 4
# for i in range(1, n):
#     for j in range(i+1, n+1):
#         for k in range(1, n+1):
#             for l in range(1, n+1):
#                 if l != k:
#                     if (dft(w_ij_kl(i, j, k, l), n)[2] == np.array([[0, 0], [0, 2]])).all():
#                         print(i, j, k, l)
                    # for mat in dft(w_ij_kl(i, j, k, l), n):
                    #     print(mat)

for m in dft_natural(w_ij(2, 2), 4):
    print(m)
                    
# for n in range(3, 10):
#     print(representation(cycle_indicator, w_ij_mat, n, n, False))
#     print(representation(excedances, w_ij_mat, n, n, False))

# print(dft(cycle_indicator, 4))
# print(dft(excedances, 4))
# for n in range(3, 5):
#     print(nullspace(T_w_ij(n)))

def mean_var(function, n: int) -> tuple:
    # Assumes function is one-local
    mats = dft(function, n)
    return mats[0][0][0]/factorial(n), np.trace(mats[1] @ np.transpose(mats[1]))

# for n in range(3, 8):
#     print(mean_var(length, n))
