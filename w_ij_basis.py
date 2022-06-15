import numpy as np
from misc import function_sum_to_vector
from misc import function_to_vector
from perm_stats import w_ij
from perm_stats import w_ij_kl
from perm_stats import inc_seq_k
from math import factorial
from permutation import Permutation
from random_walk import representation
from random_walk import w_ij_mat
from misc import vector_to_function
from misc import mats_into_vec
from dft import one_local_dft_natural

n = 4
m = np.zeros(((n-1)**2 + 1 + ((n*(n - 3))//2)**2 + (((n - 1)*(n - 2))//2)**2, factorial(n)))
count = 0
for sigma in Permutation.group(n):
    if sigma != Permutation(4, 3, 2, 1):
        funcs = []
        for i in range(1, n+1):
            for j in range(1, n+1):
                if i != j:
                    for k in range(1, n+1):
                        for  l in range(1, n+1):
                            if k != l:
                                if w_ij_kl(i, j, k, l)(sigma, n) == 1:
                                    funcs.append(w_ij_kl(i,j, k, l))
        m[count] = function_sum_to_vector(n, *funcs)
        count += 1
print(np.linalg.matrix_rank(m))


def generate_dft_m_w_ij(n):
    count = 0
    dft_m = np.zeros(((n-1)**2 + 1, (n-1)**2 + 1))
    for sigma in Permutation.group(n):
        if inc_seq_k(sigma, n-1, n) > 0:
            w_ijs = []
            for j in range(1, n+1):
                for i in range(1, n+1):
                    if w_ij(i, j)(sigma, n) == 1:
                        w_ijs.append((i, j, 1))
            dft_m[:, count] = mats_into_vec((n-1)**2 + 1, one_local_dft_natural(w_ijs, n))
            count += 1
    return dft_m

def generate_dft_m_w_ij_kl(n):
    count = 0
    dft_m = np.zeros(((n-1)**2 + 1 + ((n*(n - 3))//2)**2 + (((n - 1)*(n - 2))//2)**2, (n-1)**2 + 1 + ((n*(n - 3))//2)**2 + (((n - 1)*(n - 2))//2)**2))
    for sigma in Permutation.group(n):
        if inc_seq_k(sigma, n-2, n) > 0:
            w_ijs = []
            for j in range(1, n+1):
                for i in range(1, n+1):
                    if w_ij(i, j)(sigma, n) == 1:
                        w_ijs.append((i, j, 1))
            dft_m[:, count] = mats_into_vec((n-1)**2 + 1, one_local_dft_natural(w_ijs, n))
            count += 1
    return dft_m

print(one_local_dft_natural([(1,1,1)], n))

n = 5
dft_m = generate_dft_m_w_ij(n)

print(dft_m)

print(np.round(np.linalg.inv(dft_m) @ mats_into_vec((n-1)**2 + 1, one_local_dft_natural([(1,1, 1)], n)), 10))

# print(m)
# print(np.linalg.matrix_rank(m))


# Permutations in S4 w increasing subsequences of length 3

