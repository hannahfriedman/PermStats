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

def one_local_dft_natural(wij_indices: list, n: int) -> list:
    ''' Given a list tuples of the form (i, j, c), this function returns the first two representations 
    using Young's natural representation of the sum of cwij '''
    result = [np.zeros((1,1)), np.zeros((n-1, n-1))]
    for i, j, c in wij_indices:
        result[0] = result[0] + np.array([[math.factorial(n-1)*c]])
        result[1] = result[1] + wij_dft_n_minus_one_one(i, j, n)
    return result
        

def wij_dft_n_minus_one_one(i: int, j: int, n: int):
    ''' Return the (n-1, 1) representation of wij using YNR '''
    result = np.zeros((n-1, n-1))
    if j == 1:
        for row in range(n-1):
            val_to_add = math.factorial(n-2)
            if row == i - 2:
                val_to_add = -math.factorial(n-1)
            for col in range(n-1):
                result[row, col] = val_to_add
    else:
        for row in range(n-1):
            val_to_add = - math.factorial(n-2)
            if row == i - 2:
                val_to_add = math.factorial(n-1)
            result[row, j - 2] = val_to_add
    return result

def dft_matrix(f: Callable[[Permutation, int], int], n: int) -> np.array:
    ''' Returns the dft of f as a block diagonal matrix using YOR'''
    ft = dft(f, n)
    result = np.array(ft[0])
    for mat in ft[1:]:
        for _ in range(mat.shape[0]):
            result = np.block([
                [result, np.zeros((result.shape[0], mat.shape[1]))],
                [np.zeros((mat.shape[0], result.shape[1])), mat]
                ])
    return result

# def dft(f: Callable[[Permutation, int], int], n: int) -> list:
#     """
#     Returns the DFT of Sn with f(sigma) as the coefficient of sigma for all sigma in Sn
#     """
#     nfac = math.factorial(n)
#     sn = Permutation.group(n)
#     result = []
#     rho = matrix_rep(n)         # map from Sn to Fourier space
#     num_mats = len(rho[Permutation()])
#     for sigma in sn:
#         perm_stat = f(sigma, n)
#         if sigma == Permutation():    # for the first permutation, append matrices to list
#             for mat in rho[sigma]:
#                 result.append(perm_stat*mat)
#         else:                         # for the following permutations, add their representations to existing matrices
#             for index in range(num_mats): 
#                 result[index] = np.add(result[index], perm_stat*rho[sigma][index])
#     return adjust_zeros(result)

def dft(f: Callable[[Permutation, int], int], n: int) -> list:
    """
    Returns the DFT of Sn with f(sigma) as the coefficient of sigma for all sigma in Sn
    """
    nfac = math.factorial(n)
    sn = Permutation.group(n)
    rho = matrix_rep(n)         # map from Sn to Fourier space
    num_mats = len(rho[Permutation()])
    result = [f(Permutation(), n) * mat for mat in rho[Permutation()]]
    for sigma in sn:
        perm_stat = f(sigma, n)
        if sigma != Permutation(): # for the following permutations, add their representations to existing matrices
            for index in range(num_mats): 
                result[index] = np.add(result[index], perm_stat*rho[sigma][index])
    return adjust_zeros(result)


def inverse_dft(mats: list, n: int) -> list:
    ''' Returns the inverse dft given a list of matrices'''
    n_fac = math.factorial(n)
    reps = matrix_rep(n)
    result = []
    for g in Permutation.group(n):
        entry = 0
        for i in range(len(mats)):
            entry += mats[i].shape[0] * np.trace(reps[g.inverse()][i] @ mats[i])/n_fac
        result.append(round(entry, 10))
    return result

def inverse_dft_natural(mats: list, n: int) -> list:
    n_fac = math.factorial(n)
    reps = matrix_rep_natural(n)
    result = []
    for g in Permutation.group(n):
        entry = 0
        for i in range(len(mats)):
            entry += mats[i].shape[0] * np.trace(reps[g.inverse()][i] @ mats[i])/n_fac
        result.append(round(entry, 10))
    return result


def projection(function: Callable[[Permutation, int], int], n: int, *args) -> list:
    ''' Project f into the frequency spaces listed in args'''
    frequency = dft(function, n)
    for i in range(len(frequency)):
        if i not in args: # Replace the spaces not wanted with zeros
            frequency[i] = np.zeros((frequency[i].shape[0], frequency[i].shape[1]))
    return inverse_dft(frequency, n)

# All dft natural functions are still buggy
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
    returns a dict every permutation in Sn to it's natural matrix representation
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
 #    print(tableaux_by_shape[2])
    rho = {}
    # Loop over numbers 1 through n-1 to represent adjacent transpositions (1,2) to (n-1,n)
    for i in range(1,n):
        representation = []
        for shape in tableaux_by_shape:
            Tableau.sort(shape)
            rep = np.zeros((len(shape), len(shape))) # Representations are indexed by standard tableaux of a given shape
            for index in range(len(shape)):          # Loop over tableaux to fill in matrices
                # if i == 1 and len(shape[0].data) == 3:
                #     print(shape)
                tableau = Tableau(shape[index].data) # Create a temporary variable name for the current tableau
                col = 0
                in_same_row = False
                in_same_column = False
                for row in range(len(tableau.data)):
                    if i in tableau.data[row]:
                        col = tableau.data[row].index(i)
                        if i+1 in tableau.data[row]:
                            in_same_row = True
                            tableau_with_descent = Tableau(copy.deepcopy(tableau.data)).apply_permutation(Permutation.cycle(i, i+1))
                            for new_tableau, value in find_garnir_element(tableau_with_descent, i, i+1, row, col, tableau.data[row].index(i+1)):
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
#                    if tableau == Tableau([[1, 3, 4], [2, 5]]) and i == 2:
#                        print(tableau, tableau.switch(i), switched_index, shape)
                    rep[switched_index, index] = 1
            representation.append(rep) # Add this representation to the list of representations
        # print(i, i+1)
        # for rep in representation:
        #     print(rep)
        rho[Permutation.cycle(i, i+1)] = representation # Add representation to dictionary
    return rho

# def matrix_rep_gen_natural(n: int) -> dict:
#     """
#     n--int n in S_n
#     returns a dict that maps the generators of S_n to their orthogonal matrix representations
#     """
#     partitions = sort_partitions(generate_partitions(n))
#     tableaux_by_shape = [Tableau.gen_by_shape(n, partition) for partition in partitions]
#     rho = {}
#     # Loop over numbers 1 through n-1 to represent adjacent transpositions (1,2) to (n-1,n)
#     for i in range(1,n):
#         representation = []
#         for shape in tableaux_by_shape:
#             Tableau.sort(shape)
#             rep = np.zeros((len(shape), len(shape))) # Representations are indexed by standard tableaux of a given shape
#             for index in range(len(shape)):          # Loop over tableaux to fill in matrices
#                 tableau = Tableau(shape[index].data) # Create a temporary variable name for the current tableau
#                 col = 0
#                 in_same_row = False
#                 in_same_column = False
#                 for row in range(len(tableau.data)):
#                     if i in tableau.data[row]:
#                         col = tableau.data[row].index(i)
#                         if i+1 in tableau.data[row]:
#                             in_same_row = True
#                             col_i_plus_one = tableau.data[row].index(i+1)
#                             tableau_with_descent = Tableau(copy.deepcopy(tableau.data))
#                             tableau_with_descent.data[row][col] = i+1
#                             tableau_with_descent.data[row][col_i_plus_one] = i
#                             indices = same_row_coeff(tableau_with_descent, index, shape, i, i+1, row, col_i_plus_one, col)
#                         else:
#                             for r in range(row, len(tableau.data)):
#                                 if len(tableau.data[r]) > col and tableau.data[r][col] == i+1:
#                                     indices = same_col_coeff(tableau, index, shape, i, i+1, row, r, col)
#                                     in_same_column = True
#                             if not in_same_row and not in_same_column:
#                                 row_i_plus_one, col_i_plus_one = tableau.find(i+1)
#                                 switched = tableau.switch(i)
#                                 indices = no_sharing_coeff(switched, index, shape, i, i+1, row_i_plus_one, col_i_plus_one, row, col)
#                 for row, col, val in indices:
#                     rep[row][col] = val
#             representation.append(rep) # Add this representation to the list of representations
#         rho[Permutation.cycle(i, i+1)] = representation # Add representation to dictionary
#     return rho


# def same_row_coeff(tableau, index, shape, i, j, row, col_i, col_j):
#     a = {tableau.data[row][col_i] for row in range(row, len(tableau.data)) if len(tableau.data[row]) > col_i}
#     b = {tableau.data[row][col_j] for row in range(0, row + 1)}
#     # if len(a) == 1 and len(b) == 1:
#     #     switched_tableau = Tableau(copy.deepcopy(tableau.data))
#     #     switched_tableau.data[row_i][col_i] = j
#     #     switched_tableau.data[row_j][col_j] = i
#     #     for j in range(len(shape)):
#     #         if shape[j] == switched_tableau:
#     #             return [(index, index, - 1)]
#     c = a.union(b)
#     result = []
#     for a_prime in combinations(c, len(a)):
#         b_prime_list = sorted(list(c.difference(set(a_prime))))
#         a_prime_list = sorted(list(a_prime))
#         print(tableau)
#         print(a_prime_list, b_prime_list)
#         new_data = copy.deepcopy(tableau.data)
#         for r in range(len(new_data)):
#             if r <= row:
#                 new_data[r][col_j] = b_prime_list[r]
#             if r >= row:
#                 if len(new_data[r]) <= col_i:
#                     break
#                 new_data[r][col_i] = a_prime_list[r - row]
#         new_tableau = Tableau(new_data)
#         if new_tableau != tableau:
#             row_desc = new_tableau.row_descent()
#             print(new_tableau)
#             perm_diff_sgn = tableau.permutation_difference(new_tableau).sign
#             if row_desc != []:
#                 result += same_row_coeff(new_tableau, index, shape, *row_desc[0])
#             else:
#                 tab_no_col_desc = new_tableau.eliminate_column_descents()
#                 col_desc_factor = new_tableau.permutation_difference(tab_no_col_desc).sign
#                 for j in range(len(shape)):
#                     if shape[j] == tab_no_col_desc:
#                         result.append((j, index,  - perm_diff_sgn * col_desc_factor))
#                         break
#     return result

# def same_col_coeff(tableau, index, shape, i, j, row_i, row_j, col):
#     return [(index, index, -1)]

# def no_sharing_coeff(tableau, index, shape, i, j, row_i, row_j, col_i, col_j):
#     for i in range(len(shape)):
#         if shape[i] == tableau:
#             return [(i, index, 1)]

def find_garnir_element(tableau, i, j, row_i, col_i, col_j):
    a = {tableau.data[row][col_i] for row in range(row_i, len(tableau.data)) if len(tableau.data[row]) > col_i}
    b = {tableau.data[row][col_j] for row in range(0, row_i + 1)}
    if len(a) == 1 and len(b) ==1:
        tableau_copy = Tableau(copy.deepcopy(tableau.data))
        return [(tableau_copy.apply_permutation(Permutation.cycle(i, j)), 1)]
    c = a.union(b)
    tableaux_with_sign = []
    for a_prime in combinations(c, len(a)):
        b_prime = c.difference(set(a_prime))
#         if min(a_prime) < max(b_prime):
            # if tableau == Tableau([[1,2], [4, 3], [5]]):
            #     print(a_prime, b_prime)
        a_prime_list = sorted(list(a_prime))
        b_prime_list = sorted(list(b_prime))
        new_data = copy.deepcopy(tableau.data)
        for row in range(len(new_data)):
            if row <= row_i:
                new_data[row][col_j] = b_prime_list[row]
            if row >= row_i:
                if len(new_data[row]) <= col_i:
                    break
                new_data[row][col_i] = a_prime_list[row - row_i]
        new_tableau = Tableau(new_data)
        if new_tableau != tableau:
            row_desc = new_tableau.row_descent()
            perm_diff_sgn = tableau.permutation_difference(new_tableau).sign
            if row_desc != []:
                tableaux_with_sign += [(tab, - perm_diff_sgn * sign) for tab, sign in find_garnir_element(new_tableau, *row_desc[0])]
            else:
                tab_no_col_desc = new_tableau.eliminate_column_descents()
                col_desc_factor = new_tableau.permutation_difference(tab_no_col_desc).sign
                tableaux_with_sign.append((new_tableau, - perm_diff_sgn * col_desc_factor))
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
