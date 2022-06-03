import numpy as np
from numpy.linalg import matrix_power
import math
import copy
from itertools import combinations

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

def dft_natural(f: Callable[[Permutation, int], int], n: int) -> list:
    """
    Returns the DFT of Sn with f(sigma) as the coefficient of sigma for all sigma in Sn
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    result = []
    rho = matrix_rep_natural(n)         # map from Sn to Fourier space
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
    rho_gen = matrix_rep_gen(n)
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

def matrix_rep_natural(n: int) -> dict:
    """
    n--int n in S_n
    returns a dict every permutation in Sn to it's orthogonal matrix representation
    """
    # Create dictionary representation the DFT for the generates of Sn (adjacent transpositions)
    rho_gen = matrix_rep_gen_natural(n)
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

def matrix_rep_gen_natural(n: int) -> dict:
    """
    n--int n in S_n
    returns a dict that maps the generators of S_n to their orthogonal matrix representations
    """
    partitions = sort_partitions(generate_partitions(n))
    tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
    print(tableaux_by_shape)
    rho = {}
    # Loop over numbers 1 through n-1 to represent adjacent transpositions (1,2) to (n-1,n)
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            Tableau.sort(shape)
            rep = np.zeros((len(shape), len(shape))) # Representations are indexed by standard tableaux of a given shape
            for index in range(len(shape)):          # Loop over tableaux to fill in matrices
                tableau = Tableau(shape[index].data) # Create a temporary variable name for the current tableau
                col = 0
                in_same_row = False
                in_same_column = False
                for row in range(len(tableau.data)):
                    if i in tableau.data[row]:
                        col = tableau.data[row].index(i)
                        if i+1 in tableau.data[row]:
                            in_same_row = True
                            for new_tableau, value in find_garnir_element(tableau, i, row, col, tableau.data[row].index(i+1)):
                                new_tableau_index = 0
                                for j in range(len(shape)):     # Find the index of the resulting standard tableau
                                    if shape[j] == new_tableau:  
                                        new_tableau_index = j
                                        break
                                rep[new_tableau_index, index] = value
                        break
                if not in_same_row:
                    for row in tableau.data:
                        if len(row) > col and row[col] == i+1:
                            rep[index, index] = -1
                            in_same_column = True
                if not in_same_row and not in_same_column:
                    switched_index = 0
                    for j in range(len(shape)):     # Find the index of the resulting standard tableau
                        if shape[j] ==  tableau.switch(i):  
                            switched_index = j
                            break
                    rep[switched_index, index] = 1
            representation.append(rep) # Add this representation to the list of representations
        print(i, i+1)
        for rep in representation:
            print(rep)
        rho[Permutation.cycle(i, i+1)] = representation # Add representation to dictionary
    return rho

def find_garnir_element(tableau, i, row_i, col_i, col_i_plus_one):
    a = {tableau.data[row][col_i] for row in range(row_i, len(tableau.data)) if len(tableau.data[row]) > col_i}
    b = {tableau.data[row][col_i_plus_one] for row in range(0, row_i + 1)}
    if len(a) == 1 and len(b) ==1:
        return [(tableau, 1)]
    c = a.union(b)
    tableaux_with_sign = []
    for a_prime in combinations(c, len(a)):
        b_prime = c.difference(set(a_prime))
        if min(a_prime) < max(b_prime):
            a_prime_list = sorted(list(a_prime))
            b_prime_list = sorted(list(b_prime))
            new_data = copy.deepcopy(tableau.data)
            for row in range(len(new_data)):
                if row <= row_i:
                    new_data[row][col_i_plus_one] = b_prime_list[row]
                if row >= row_i:
                    if len(new_data[row]) <= col_i:
                        break
                    new_data[row][col_i] = a_prime_list[row - row_i]
            new_tableau = Tableau(new_data)
            tableaux_with_sign.append((new_tableau, tableau.permutation_difference(new_tableau).sign))
    print(i)
    for pair in tableaux_with_sign:
        print(pair)
    return tableaux_with_sign

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
