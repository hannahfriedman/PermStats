import numpy as np
from numpy.linalg import matrix_power
import math
import copy

from tableau import Tableau
from permutation import Permutation
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

def matrix_rep(n: int) -> dict:
    """
    n--int n in S_n
    returns a dict every permutation in Sn to it's orthogonal matrix representation
    """
    # Create dictionary representation the DFT for the generates of Sn (adjacent transpositions)
    rho_gen =(n)
    rho = {}
    # Creates a list of lists of permutation in Sn factored into adjacent transpositions
    sn = factor_sn(n)
    for sigma in sn:
        key = Permutation()
        val = []
        # Becuase the identity is not in the list of generators, start by populating val list 
        # with (1 2)^2 which is the identity
        for mat in rho_gen[Permutation.cycle(1,2)]:
            val.append(matrix_power(mat, 2))
        # For each factor, multiply by the appropriate matrix
        for transposition in sigma:
            if transposition != Permutation(): # don't use the identity because we don't have a matrix representation for it
                key *= transposition
                for i in range(len(val)):
                    val[i] = val[i] @ rho_gen[transposition][i]
        rho[key] = val # Add key and value to map
    return rho


def factor_sn(n: int) -> list:
    """
    n - int n in S_n
    returns a list of lists of factorizations of all elements of Sn
    We find this list recursively but multiplying Sn-1 by transpositions that move n into
    all possible positions
    """
    # Base case: permutations in Sn are already factored
    if n == 2:
        return [[Permutation.cycle()], [Permutation.cycle(1,2)]]
    else:
        # Find factorization of Sn-1 recursively
        sn_minus_one = factor_sn(n-1)
        # Sn-1 is a subset of Sn, so start by adding those values to our list
        sn = sn_minus_one
        prev_coset = sn_minus_one
        curr_coset = []
        for i in range(1, n):
            # Move n to next position in permutation for each sigma in the previous coset
            for sigma in prev_coset:
                tau = sigma + [Permutation.cycle(n-i,n-i+1)]
                curr_coset.append(tau)
            # Add the values we just computed to sn
            sn = sn + curr_coset
            # Reset necessary values for next iteration
            prev_coset = curr_coset
            curr_coset = []
        return sn        

def matrix_rep_gen(n: int) -> dict:
    """
    n--int n in S_n
    returns a dict that maps the generators of S_n to their orthogonal matrix representations
    """
    partitions = sort_partitions(generate_partitions(n))
    tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
    rho = {}
    # Loop over numbers 1 through n-1 to represent adjacent transpositions (1,2) to (n-1,n)
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            Tableau.sort(shape)
            rep = np.zeros((len(shape), len(shape))) # Representations are indexed by standard tableaux of a given shape
            for index in range(len(shape)):          # Loop over tableaux to fill in matrices
                tableau = Tableau(shape[index].data) # Create a temporary variable name for the current tableau
                rep[index, index] = 1/(tableau.signed_distance(i))  # Diagonal matrix entries
                switched = tableau.switch(i)         # Apply the transposition to the tableau
                if switched.is_standard():           # If the result is a standard tableau...
                    switched_index = 0
                    for j in range(len(shape)):     # Find the index of the resulting standard tableau
                        if shape[j] == switched:  
                            switched_index = j
                            break
                    rep[switched_index, index] = math.sqrt(1 - (shape[index].signed_distance(i))**(-2)) # Fill the entry indexed by the current tableau and the switched tableau
            representation.append(rep) # Add this representation to the list of representations
        rho[Permutation.cycle(i, i+1)] = representation # Add representation to dictionary
    return rho

def generate_partitions(n: int) -> list:
    '''
    Generates all partitions of size n
    Returns a list of lists, showing the size of each partition
    Works recursively
    '''
    ans = []
    # Base case 1
    if n == 1:
        return [[1]]
    # Base case 2
    elif n == 0:
        return [[]]
    # Convert all partitions of integers smaller than n to partitions of n by adding the 
    # necessary value
    for x in range(1, n):
        ans += [[x] + part for part in generate_partitions(n-x)]
    # Remove duplicates and append the partition [n], which is not accounted for above
    return remove_dubs(ans) + [[n]]

def sort_partitions(partitions: list) -> list:
    '''
    sort partitions according to dominance order, in descending order using bubble sort
    '''
    switched = True
    while switched:
        switched = False
        # Loop over all partitions
        for i in range(len(partitions)-1):
            # If  adjacent partitions are not in the right order, switch them
            if compare_partitions(partitions[i+1], partitions[i]): 
                partitions[i], partitions[i+1] = partitions[i+1], partitions[i]
                switched = True # We haven't finished sorting yet
    return partitions

def compare_partitions(partition1: list, partition2: list) -> bool:
    '''
    returns true if partition1 dominates partition2, false otherwise
    '''
    sum1 = 0
    sum2 = 0
    for i in range(len(partition1)):
        sum1 += partition1[i]
        sum2 += partition2[i]
        if sum1 != sum2:
            return sum1 > sum2

def remove_dubs(partition: list) -> list:
    ''' 
    Removes duplicates in a list of lists
    Makes sure that any inner lists are sorted first (treats inner lists as multi-sets)
    ''' 
    # Sort all the partitions and then put them in descending order
    for part in partition:
        part.sort()
        part.reverse()
    # Adds only unique partitions to result and returns
    result = []
    for part in partition:
        if part not in result:
            result.append(part)
    return result
