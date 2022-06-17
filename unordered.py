from perm_stats import w_ij_kl_unordered
import numpy as np
from dft import dft
from random_walk import representation
from random_walk import w_ij_kl_unordered_mat
from random_walk import w_ij_kl_mat
from misc import choose
from permutation import Permutation
from decompress import T_w_ij
from decompress import T_w_ij_kl
from random_walk import rho_reg
from math import factorial
from perm_stats import w_ij_kl
from perm_stats import w_ij
from perm_stats import fixed_points
from perm_stats import pairs_of_fixed_points
from perm_stats import unordered_fixed_points
from misc import function_to_vector
from misc import vector_to_function
from itertools import combinations
from perm_stats import length
from perm_stats import major_index
from perm_stats import descent
from perm_stats import excedances
from dft import projection
i = 3
j = 2
k = 1
l = 4
n = 6

def indicator(sigma, n):
    return int(sigma == Permutation(4, 2, 3, 1))

# print(dft(w_ij_kl_unordered(i, j, k, l), 4))
# print(representation(indicator, w_ij_kl_unordered_mat, choose(n, 2), n, False))


# print(np.transpose(T_w_ij(n)) @ T_w_ij(n))
# print(np.transpose(T_w_ij_kl(n)) @ T_w_ij_kl(n))
n = 4
v_sum = np.zeros((factorial(n)))
reg_rep = np.zeros((factorial(n), factorial(n)))
for sigma in Permutation.group(n):
    for i in range(1, n+1):
        for j in range(1, n+1):
            if i != j:
                reg_rep += w_ij_kl(i, j, i, j)(sigma, n) * rho_reg(sigma, n)
                v_sum = v_sum + (w_ij_kl(i, j, i, j)(sigma, n) * np.array(function_to_vector(w_ij_kl(i, j, i, j), n)))
# print(dft(vector_to_function(v_sum, n), n))
    
print(np.linalg.eig(reg_rep)[0])
# print((reg_rep == np.transpose(T_w_ij_kl(n)) @ T_w_ij_kl(n)).all())
# print(representation(pairs_of_fixed_points, w_ij_kl_mat, n * (n-1), n, False))

for n in range(4, 6):
    # print(dft(fixed_points, n)[1])
    print(dft(unordered_fixed_points, n)[1])
    print(dft(pairs_of_fixed_points, n)[0])
    print(dft(pairs_of_fixed_points, n)[1])
    print(dft(pairs_of_fixed_points, n)[2])
    print(dft(pairs_of_fixed_points, n)[3])
    

def T_w_ij_kl_unordered(n):
    result = np.zeros((choose(n, 2)**2, factorial(n)))
    row = 0
    for i, j in combinations(list(range(1, n+1)), 2):
        for k, l in combinations(list(range(1, n+1)), 2):
            result[row] = function_to_vector(w_ij_kl_unordered(i, j, k, l), n)
            row += 1
    return result
            
            
def projection_one(f, n):
    t1t1star = np.ones((factorial(n), factorial(n)))
    return t1t1star @ f /factorial(n)


def projection_two(f, n):
    t1t1star = np.ones((factorial(n), factorial(n)))
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    return (t2t2star - t1t1star) @ f * (n - 1)/factorial(n)


def projection_three(f, n):
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    t3t3star = np.transpose(T_w_ij_kl_unordered(n)) @ T_w_ij_kl_unordered(n)
    return (t3t3star - t2t2star) @ f * (n - 3)/(2 * factorial(n - 1))


def projection_four(f, n):
    t1t1star = np.ones((factorial(n), factorial(n)))
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    t3t3star = np.transpose(T_w_ij_kl_unordered(n)) @ T_w_ij_kl_unordered(n)
    t4t4star = np.transpose(T_w_ij_kl(n)) @ T_w_ij_kl(n)
#     print(t4t4star - t3t3star/2 - t2t2star/2 + t1t1star/2)
    return (t4t4star - t3t3star - t2t2star + t1t1star) @ f/ (2 * n * factorial(n - 3))
    # return (t4t4star - t3t3star/2 - t2t2star/2 + t1t1star/2)@ f/(n * factorial(n-3))
#     return (t4t4star @ f - projection_one(f, n) * factorial(n) / 2 - projection_two(f, n) * factorial(n)/(n-1) - projection_three(f, n) * factorial(n-1)/(n - 3))/(n * factorial(n - 3))

n = 4
print(np.linalg.eig(np.transpose(T_w_ij_kl(n)) @ T_w_ij_kl(n))[0])

f = length
v = function_to_vector(f, n)
print(v)
print(projection_one(v, n))
print(projection_two(v, n))
print(projection_three(v, n))
print(projection_four(v, n))
# print(dft(vector_to_function(projection_four(v, n), n), n))
print(projection(length, n, 3))
# print(dft(pairs_of_fixed_points, n))
triv = T_w_ij_kl(n) @ projection_one(v, n)/factorial(n)
n_min_one_one = T_w_ij_kl(n) @ projection_two(v, n)*(n-1)/(factorial(n) * 2)
n_min_two_two = T_w_ij_kl(n) @ projection_three(v, n) * (n-3)/(factorial(n-1) * 2)
n_min_two_one_one = T_w_ij_kl(n) @ projection_four(v, n) / (2 * n * factorial(n-3))
print(v)
print(np.transpose(T_w_ij_kl(n))  @ (triv + n_min_one_one + n_min_two_two + n_min_two_one_one))
print(triv + n_min_one_one + n_min_two_two + n_min_two_one_one)
                                                              
