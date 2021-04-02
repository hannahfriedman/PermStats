import numpy as np
from numpy.linalg import matrix_power
import math
import copy
from Tableau import Tableau
from permutation import Permutation
from misc import matlab_syntax
from misc import adjust_zeros
from typing import Callable
import perm_stats

#Docs for Permutation class: 
#https://permutation.readthedocs.io/en/stable/_modules/permutation.html#Permutation.cycle


def dft(f: Callable[[Permutation, int], int], n: int) -> list:
    """
    Returns the DFT of Sn with f(sigma) as the coefficient of sigma for all sigma in Sn
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    result = []
    rho = matrix_rep(n)         # map from Sn to Fourier space
    num_mats = len(rho[Permutation()])
    for sigma in sn:
        perm_stat = f(sigma, n)
        if sigma == Permutation():    # for the first permutation, append matrices to list
            for mat in rho[sigma]:
                result.append(perm_stat*mat)
        else:                         # for the following permutations, add their representations to existing matrices
            for index in range(num_mats): 
                result[index] = np.add(result[index], perm_stat*rho[sigma][index])
    return adjust_zeros(result)

def matrix_rep(n):
    """
    n--int n in S_n
    returns a dict that maps the 2-cycles of S_n to their orthogonal matrix representations
    """
    rho_gen = matrix_rep_gen(n)
    rho = {}
    sn = factor_sn(n)
    for perm in sn:
        key = Permutation.cycle()
        val = []
        for mat in rho_gen[Permutation.cycle(1,2)]:
            val.append(matrix_power(mat, 2))
        for transposition in perm:
            if transposition != Permutation.cycle():
                key *= transposition
                for i in range(len(val)):
                    val[i] = np.matmul(val[i], rho_gen[transposition][i])
        rho[key] = val
    return rho


def factor_sn(n):
    """
    n--int n in S_n
    returns a list of lists of factorizations of all elements of Sn
    """
    if n == 2:
        return [[Permutation.cycle()], [Permutation.cycle(1,2)]]
    else:
        sn_minus_one = factor_sn(n-1)
        sn = sn_minus_one
        prev_coset = sn_minus_one
        curr_coset = []
        for i in range(1, n):
            for perm in prev_coset:
                newPerm = perm + [Permutation.cycle(n-i,n-i+1)]
                curr_coset.append(newPerm)
            sn = sn + curr_coset
            prev_coset = curr_coset
            curr_coset = []
        return sn        

def matrix_rep_gen(n):
    """
    n--int n in S_n
    returns a dict that maps the generators of S_n to their orthogonal matrix representations
    """
    partitions = generate_partitions(n)
    partitions.reverse()
    tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
    rho = {}
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            Tableau.sort(shape)
            rep = np.zeros((len(shape), len(shape)))
            for index in range(len(shape)):
                tableau = Tableau(shape[index].data)
                rep[index, index] = 1/(tableau.signed_distance(i))
                switched = tableau.switch(i)
                if switched.is_standard():
                    switched_index = 0
                    for j in range(len(shape)):
                        if shape[j] == switched:
                            switched_index = j
                            break
                    rep[switched_index, index] = math.sqrt(1 - (shape[index].signed_distance(i))**(-2))
            representation.append(rep)
        rho[Permutation.cycle(i, i+1)] = representation
    return rho

def generate_partitions(n):
    '''
    Generates all partitions of size n
    Returns a list of lists, showing the size of each partition
    '''
    ans = []
    if n == 1:
        return [[1]]
    elif n == 0:
        return [[]]
    for x in range(1, n):
        ans += [[x] + part for part in generate_partitions(n-x)]
    return remove_dubs(ans) + [[n]]

def remove_dubs(partition):
    ''' 
    Removes duplicates in a list of lists
    Makes sure that any inner lists are sorted first (treats inner lists as multi-sets)
    ''' 
    for part in partition:
        part.sort()
        part.reverse()
    result = []
    for part in partition:
        if part not in result:
            result.append(part)
    return result


def matrix_rep_transpositions(n):
    """
    n--int n in S_n
    returns a dict that maps the 2-cycles of S_n to their orthogonal matrix representations
    With help from formulas obtained in: https://math.stackexchange.com/questions/3420570/writing-permutations-as-products-of-adjacent-transposition 
    """
    rho_gen = matrix_rep_gen(n)
    rho = rho_gen
    
    # Calcultes the matrix representation for all 2-cycles
    for diff in range(2, n):
        for startVal in range(1, n - diff + 1):
            endVal = startVal + diff
            permutation = Permutation.cycle(startVal, endVal)
            perm_factor1 = Permutation.cycle(endVal - 1, endVal)
            perm_factor2 = Permutation.cycle(startVal, endVal-1)
            matrix_factor1 = rho[perm_factor1]
            matrix_factor2 = rho[perm_factor2]
            matrix_rep = [np.matmul(np.matmul(matrix_factor1[i], 
                                    matrix_factor2[i]), 
                                    matrix_factor1[i])
                                    for i in range(len(matrix_factor1))]
            rho[permutation] = matrix_rep
    
    # Matrix representation for the rest
    return rho

def __main__():
    """
    rho = matrix_rep(3)
    print(rho)
    sn = Permutation.group(3)
    for sigma in sn:
    print(rho[sigma][1])"""
    dft(perm_stats.excedances, 3)


__main__()