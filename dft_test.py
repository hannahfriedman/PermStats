from dft import dft
from dft import dft_matrix
from dft import dft_natural
from dft import inverse_dft
from dft import projection
from perm_stats import w_ij_kl
from perm_stats import w_ij
from perm_stats import cycle_indicator
from random_walk import representation
from random_walk import w_ij_mat
from perm_stats import excedances
from perm_stats import length
from perm_stats import descent
from perm_stats import major_index
from decompress import T_w_ij
import numpy as np
from misc import nullspace
from misc import function_to_vector
from math import factorial
from misc import left_action_on_vector
from misc import left_right_action_on_vector
from small_dft import n_minus_one_one
from permutation import Permutation
from misc import vector_to_function
import random
from one_local_stats import f
from one_local_stats import compute_distribution
from random_walk import rho_reg


def all_even(mat):
    for i in mat:
        for j in i:
            if j%2 == 1:
                return False
    return True


# n = 4
# ITER = 1000
# sampling_space = list(range(factorial(n)))
# num_one_local = 0
# for run in range(ITER):
#     function = factorial(n) * [1]
#     zero = random.randint(0, factorial(n) - 1)
#     three = random.randint(0, factorial(n) - 1)
#     while zero == three:
#         three =	random.randint(0, factorial(n) - 1)
#     function[three] = 3
#     function[zero] = 0
#     if tuple(function) in f():
#         print('yay')
#     num_twos = 0
#     while num_twos < 11:
#         entry = random.choice(sampling_space)
#         if entry != three and entry != zero and function[entry] != 2:
#             function[entry] = 2
#             num_twos += 1
#     ft = dft(vector_to_function(function, n), n)
#     one_local = True
#     for i in range(2, len(ft)):
#         zeros = np.zeros((ft[i].shape[0], ft[i].shape[1]))
#         if not (ft[i] == zeros).all():
#             one_local = False
#             break
#     if not one_local and all_even(representation(vector_to_function(function, n), w_ij_mat, n, n, False)):
#         print(function, representation(vector_to_function(function, n), w_ij_mat, n, n, False))
#     num_one_local += int(one_local)
# print(num_one_local/1000)


# n = 4    
# funcs = f()
# index = 500
# for i in range(len(funcs)):
#     if compute_distribution(funcs[i], 4) == compute_distribution(funcs[index], 4): # [1, 11, 11, 1, 0]:
#         print(funcs[i], '\n', representation(vector_to_function(funcs[i], n), w_ij_mat, n, n, False))
#         print(compute_distribution(funcs[index], 4))

# n = 4    
# funcs = f()
# index = 500
# for i in range(len(funcs)):
#     if compute_distribution(funcs[i], 4) == [1, 11, 11, 1, 0]:
#         print(funcs[i], '\n', representation(vector_to_function(funcs[i], n), w_ij_mat, n, n, False))
#         print(compute_distribution(funcs[index], 4))

# func = vector_to_function((1, 2, 0, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 1, 1, 2, 2, 1, 1), 4)
def temp(sigma, n):
    return func(Permutation(1, 4, 2, 3)*sigma * Permutation(1, 2, 4, 3), n)

def indicator(sigma, n):
    return int(sigma == Permutation(3, 2, 1, 4))
n = 4
# print(function_to_vector(temp, n))
# print(representation(temp, w_ij_mat, n, n, False))

transform = dft_matrix(indicator, n)
regular = rho_reg(Permutation(3, 2, 1, 4), n)
print((transform == np.linalg.inv(transform) @ regular @ transform).all())

# n = 4
# print(tuple(function_to_vector(excedances, n)) in f())



# n = 3
# vector = [2, 1, 1, 1, 1, 0]
# f = vector_to_function(vector, n)
# print(projection(f, n, 0, 1))

# i = 2
# j = 1
# k = 4
# l = 1
# n = 4
# function = length
# print(np.array(projection(function, n, 1)), np.array(function_to_vector(function, n)))

# print(projection(w_ij(n, n), n, 1))
# print(left_action_on_vector(projection(w_ij(n, n), n, 1), [g for g in Permutation.group(n)], Permutation.cycle(n, n-1)))
# print(left_action_on_vector(projection(w_ij(n, n), n, 1), [g for g in Permutation.group(n)], Permutation.cycle(n, n-2)))
# print(left_action_on_vector(projection(w_ij(n, n), n, 1), [g for g in Permutation.group(n)], Permutation.cycle(n, n-3)))
# print(left_right_action_on_vector(projection(w_ij(n, n), n, 1), [g for g in Permutation.group(n)], Permutation.cycle(n, n-3), Permutation.cycle(n, n-3)))
#print(dft_natural(w_ij_kl(i, j, k, l), n)[1])
# for i in range(1, n):
#     for j in range(i+1, n+1):
#         for k in range(1, n+1):
#             for l in range(1, n+1):
#                 if l != k:
#                     print(i, j, k, l)
#                     print(dft_natural(w_ij_kl(i, j, k, l), n)[3])
                    # for mat in dft_natural(w_ij_kl(i, j, k, l), n):
                    #     print(mat)

# for m in dft_natural(w_ij_kl(2, 4, 1, 2), 4):
#     print(m)
                    
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
    return mats[0][0][0]/factorial(n), np.trace(mats[1] @ np.transpose(mats[1])/ factorial(n-1))

# for n in range(3, 8):
#     print(mean_var(length, n))
# for n in range(3, 8):
#     print(mean_var(excedances, n))
