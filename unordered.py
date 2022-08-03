from perm_stats import w_ij_kl_unordered
import numpy as np
from dft import dft
from random_walk import representation
from random_walk import rep
from random_walk import rep_w_ij
from random_walk import w_ij_kl_unordered_mat
from random_walk import w_ij_kl_mat
from random_walk import w_ij_mat
from misc import choose
from permutation import Permutation
from decompress import T_w_ij
from decompress import T_w_ij_fast
from decompress import T_w_ij_rec
from decompress import T_w_ij_rec_fast
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
from decompress import T_w_ij_kl_natural
from decompress import T_w_ij_natural
import cProfile

def T_w_ij_kl_unordered(n):
    ''' Analysis operator into M^(n - 2, 2) space '''
    result = np.zeros((choose(n, 2)**2, factorial(n)))
    row = 0
    for i, j in combinations(list(range(1, n+1)), 2):
        for k, l in combinations(list(range(1, n+1)), 2):
            result[row] = function_to_vector(w_ij_kl_unordered(i, j, k, l), n)
            row += 1
    return result


def projection_one(f, n):
    ''' Compute trivial projection '''
    t1t1star = np.ones((factorial(n), factorial(n)))
    return t1t1star @ f /factorial(n)

def projection_two(f, n):
    ''' Compute projection into S^(n - 1, 1) '''
    t1t1star = np.ones((factorial(n), factorial(n)))
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    return (t2t2star - t1t1star) @ f * (n - 1)/factorial(n)


def projection_three(f, n):
    ''' Compute projection into S^(n - 2, 2) '''
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    t3t3star = np.transpose(T_w_ij_kl_unordered(n)) @ T_w_ij_kl_unordered(n)
    return (t3t3star - t2t2star) @ f * (n - 3)/(2 * factorial(n - 1))

def projection_four(f, n):
    ''' Compute projection into S^(n - 2, 1, 1) '''    
    t1t1star = np.ones((factorial(n), factorial(n)))
    t2t2star = np.transpose(T_w_ij(n)) @ T_w_ij(n)
    t3t3star = np.transpose(T_w_ij_kl_unordered(n)) @ T_w_ij_kl_unordered(n)
    t4t4star = np.transpose(T_w_ij_kl(n)) @ T_w_ij_kl(n)
    return (t4t4star - t3t3star - t2t2star + t1t1star) @ f/ (2 * n * factorial(n - 3))

def recover_w_ij_kl_coefficients(f, n):
    ''' Given a function, find the wijkl coefficients needed to write the two-local projection of f'''
    v = function_to_vector(f, n)
    return T_w_ij_kl(n) @ projection_one(v, n)/factorial(n) + T_w_ij_kl(n) @ projection_two(v, n)*(n-1)/(factorial(n) * 2) + T_w_ij_kl(n) @ projection_three(v, n) * (n-3)/(factorial(n-1) * 2) + T_w_ij_kl(n) @ projection_four(v, n) / (2 * n * factorial(n-3))

def reecover_w_ij_coefficients(f, n):
    ''' Given a function, find the wijkl coefficients needed to write the two-local projection of f'''    
    trival_proj = len(wijs) * np.ones((factorial(n), 1)) / n
    w_ij_proj = (n-1) * (np.transpose(T_w_ij_rec(n)) @ sum([w_ij_image(*wij, n) for wij in wijs])/factorial(n) - trival_proj)
    return trival_proj + w_ij_proj

def recover_one_local_function(wijs, n):
    ''' Returns the one-local projection of a functin (hopefully quicker than a dft)'''
    trival_proj = len(wijs) * np.ones((factorial(n), 1)) / n
    w_ij_proj = (n-1) * (np.transpose(T_w_ij_rec(n)) @ sum([w_ij_image(*wij, n) for wij in wijs])/factorial(n) - trival_proj)
    return trival_proj + w_ij_proj


def w_ij_image(i, j, n):
    ''' Returns a vectorized version of the permutation representation of wij '''
    cross = np.zeros((n, n))
    for index in range(n):
        cross[i-1, index] = - factorial(n-2)
        cross[index, j - 1] = - factorial(n-2)
    cross[i-1, j-1] += factorial(n - 1) + factorial(n-2)
    return np.reshape(factorial(n-2) * np.ones((n, n)) + cross, (n**2, 1))
